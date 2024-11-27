function [dist, signals] = RatioDist(A0, Ab, sigma0, opts)

% Function to construct ratio distribution of two Rician distributions
% (different means but same sigmas)

arguments

    A0 % b=0 mean signal
    Ab % b\=0 mean signal
    sigma0 % Noise standard deviation at TE=0

    % options
    opts.N0 = 1 % NSA for b=0 image
    opts.Nb = 1 % NSA for b\=0 image


    % Range/Resolution of output pdf
    opts.zs 
    opts.zmin = 0
    opts.zmax = 1
    opts.dz = 0.005

    % Range/Resolution of integral
    opts.ys
    opts.ymin = 0 
    opts.ymax = 1
    opts.dy = 0.005 

    % Calibration
    opts.calibrationfolder = fullfile('Noise Models', 'Ratio', 'Calibration Curves')
end


% == Construct array of zs

if ~isfield(opts,'zs')
    zs = linspace(opts.zmin, opts.zmax, ceil( (opts.zmax-opts.zmin)/opts.dz) +1);
else
    zs = opts.zs;
end

% == Construct array of ys (to integrate over at each z)

if ~isfield(opts,'ys')
    ys = linspace(opts.ymin, opts.ymax, ceil( (opts.ymax-opts.ymin)/opts.dy)+1);
else
    ys = opts.ys;
end


%% Construct first Rice distribution (b=0)

b0sigma = (1/sqrt(opts.N0))*sigma0;

RiceDist0 = makedist('Rician','s',A0,'sigma',b0sigma);



%% Construct second Rice distribution (b>0)

% Calibrate mean and sigma after extra averaging (low SNR)
bSNR = Ab/sigma0;

if and(bSNR<=5, opts.Nb>1)

    load(fullfile(char(opts.calibrationfolder) , ['NSA ' num2str(opts.Nb)] , 'MeanCalibration.mat'))
    SNRs = keys(MeanDict);
    bSNRrounded = roundtowardvec(bSNR, SNRs);
    Ab = Ab*MeanDict(bSNRrounded);

    load(fullfile(char(opts.calibrationfolder) , ['NSA ' num2str(opts.Nb)] , 'SigmaCalibration.mat'))
    bsigma = sigma0*SigmaDict(bSNRrounded);

else
    bsigma = (1/sqrt(opts.Nb))*sigma0;
end

RiceDistb = makedist('Rician','s',Ab,'sigma',bsigma);


%% Make Ratio distribution

% Initialise array for ratio distribution
dist = zeros(length(zs),1);

% For each z value, integrate over y array
for zindx = 1:length(zs)

    z = zs(zindx);

    % Construct integrand
    integrand = (ys).*(RiceDistb.pdf(z*ys)).*(RiceDist0.pdf(ys));

    % Integrate integrad over ys
    integral = trapz(opts.dy,integrand);

    dist(zindx) = integral;

end

signals = zs;

end