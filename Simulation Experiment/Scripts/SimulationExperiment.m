% MATLAB code to set up and run simulation experiment

% Define rootfolder
rootfolder = pwd;


%% Simulation settings

% VERDICT model
modeltype = 'Original VERDICT';

% DW-MRI imaging protocol
schemename = 'Original';

% Sigma0 value (run each sigma0 value as separate experiment)
sigma0 = 0.04;

% Fitting type
fittingtype = 'normal';  % 'normal' using Rice network; 'adaptive' using Ratio networks

if strcmp(fittingtype, 'normal')
    fittingtechnique = 'MLP';
    noisetype = 'Rice';
    T2train = 10000;
    sigma0train = 0.04;
else
    noisetype = 'Ratio';
end


%% Define grid of voxels 

% Number of voxels per grid section (side length)
ngrid = 20;

% fIC values
fICmin = 0.05;
fICmax = 0.95;
NfIC = 19;
fICs = linspace(fICmin, fICmax, NfIC);

% T2 values
T2s = [50, 60, 80];
NT2 = length(T2s);

% sigma0 values
sigma0s = sigma0*ones(size(T2s)); 
Nsigma0 = length(sigma0s);

% fIC grid
fICgrid = zeros(NT2*ngrid, NfIC*ngrid);
for fICindx = 1:NfIC
    for gridindx = 1:ngrid
        fICgrid( :, (fICindx-1)*ngrid + gridindx) = fICs(fICindx);
    end
end


% T2 grid
T2grid = zeros(NT2*ngrid, NfIC*ngrid);
for T2indx = 1:NT2
    for gridindx = 1:ngrid
        T2grid( (T2indx-1)*ngrid + gridindx, :) = T2s(T2indx);
    end
end


% sigma0 grid
sigma0grid = zeros(Nsigma0*ngrid, NfIC*ngrid);
for sigma0indx = 1:Nsigma0
    for gridindx = 1:ngrid
        sigma0grid( (sigma0indx-1)*ngrid + gridindx, :) = sigma0s(sigma0indx);
    end
end

% Intracellular comparment radii 
muRs = [7 8];
sigmaRs = [1 1];
Rs=linspace(0.1,15.1,16);


%% Simulate signals over grid

schemesfolder = fullfile(rootfolder, 'Simulation Experiment', 'Schemes');

load(fullfile(schemesfolder, [schemename '.mat']))
nscheme = length(scheme);

% Signal grid
signalgrid = zeros(NT2*ngrid, NfIC*ngrid, nscheme);
normsignalgrid = zeros(NT2*ngrid, NfIC*ngrid, nscheme);

% Simulate signals over grid
for jndx = 1:size(T2grid,1)
    for indx = 1:size(fICgrid, 2)

        disp(['Grid indices: ' num2str(jndx) ', ' num2str(indx) ])
        
        sigma0 = sigma0grid(jndx, indx);
        T2 = T2grid(jndx, indx);

        % Rs
        muR = muRs(1) + (muRs(2)-muRs(1))*rand();
        sigmaR = sigmaRs(1) + (sigmaRs(2)-sigmaRs(1))*rand();
        fRs = normpdf(Rs, muR, sigmaR);
        fRs = (1/sum(fRs))*fRs;

        % fIC
        fIC = fICgrid(jndx, indx);
        
        % fVASC
        fVASC = min(0.1, 1-fIC)*rand();

        % fEES
        fEES = 1-fIC-fVASC;


        for schemeIndx = 1:length(scheme)

            scan_params = [scheme(schemeIndx).delta, scheme(schemeIndx).DELTA, scheme(schemeIndx).bval];

            if scan_params(3)==0
                normsignalgrid(jndx,indx, schemeIndx) = 1;
                continue;
            end

            TE = scheme(schemeIndx).TE;
            NSA = scheme(schemeIndx).NSA;
            Rav = scheme(schemeIndx).Rav;

            % Tissue parameter vector
            tps = [fIC*fRs, fEES ,fVASC];
            
            % Diffusion weighting fraction
            fd = simulateSignal( ...
                tps, ...
                scan_params, ...
                modeltype,...
                Rs = Rs...
                );


            % == b=0 signal distribution

            % Noiseless signal
            b0signal = exp(-TE/T2);

            % Define Rice distribution
            b0dist = makedist('Rician', s=b0signal, sigma = sigma0 );
            
            signals = linspace(0,2,201);
            b0pdf = b0dist.pdf(signals);


            % == b>0 signal distribution

            % Noiseless signal
            bsignal = b0signal*fd;

            % Define Rice distribution
            bdist = makedist('Rician', s=bsignal, sigma = sigma0 );
            
            signals = linspace(0,2,201);
            bpdf = bdist.pdf(signals);


            % == b=0 signal

            % Take NSA samples
            sample = 0;
            for s = 1:NSA
                sample = sample + sampleDistribution(b0pdf, signals);
            end
            b0sample = sample/NSA;  

            signalgrid(jndx, indx, schemeIndx-1) = b0sample;


            % == b>0 signal

            % Take Nb samples
            Nb = NSA*Rav;
            sample = 0;
            for s = 1:Nb
                sample = sample + sampleDistribution(bpdf, signals);
            end
            bsample = sample/Nb;    

            signalgrid(jndx, indx, schemeIndx) = bsample;

            % == Normalised signal
            normsignalgrid(jndx, indx, schemeIndx) = bsample/b0sample;

        end
    end
end




%% APPLY FITTING

% Input data array for fitting
Y = normsignalgrid;

modelfolder = fullfile(pwd, 'Simulation Experiment', 'MLP', 'models');

switch fittingtype

    case 'normal'

        thismodelfolder = fullfile(modelfolder, modeltype, schemename, noisetype, ['T2_' num2str(T2train)], ['sigma_' num2str(sigma0train)]);
        
        [fIC, fEES, fVASC, R, rmse] = verdict_fit( ...
            Y, ...
            scheme, ...
            modeltype = modeltype,...
            fittingtechnique=fittingtechnique,...    
            modelfolder=thismodelfolder...
            );


    case 'adaptive'

        % Load possible T2 and sigma0 values
        T2vals = load(fullfile(modelfolder , modeltype , schemename , noisetype, 'T2s.mat')).T2s;
        sigma0vals = load(fullfile(modelfolder , modeltype , schemename , noisetype, 'sigma0s.mat')).sigma0s;

        [fIC, fEES, fVASC, R, rmse] = adaptive_fit( ...
            Y,...
            schemename,...
            modeltype = modeltype,...
            modelfolder = modelfolder,...
            schemesfolder = schemesfolder,...
            noisetype = 'Ratio',...
            T2 = T2grid,...
            sigma0 = sigma0grid,...
            T2vals = T2vals,...
            sigma0vals = sigma0vals...
            );


end


%% Display bias and std results

% Difference maps
Diffs = fIC-fICgrid;

% Biases
Biases = zeros(NT2, NfIC);

% Variances
Variances = zeros(NT2, NfIC);

for fICindx = 1:NfIC
    for T2indx = 1:NT2

        vals = Diffs((T2indx-1)*ngrid+1:T2indx*ngrid, (fICindx-1)*ngrid+1:fICindx*ngrid);
        bias = mean(vals(:));
        vari = var(vals(:));
        Biases(T2indx, fICindx) = bias;
        Variances(T2indx, fICindx) = vari;

    end
end


% BIAS FIGURE
f1=figure;
for T2indx = 1:NT2
    T2 = T2s(T2indx);
    plot(fICs, Biases(T2indx,:), '-*', DisplayName = ['T2 = ' num2str(T2) ' ms'])
    hold on
end
ylim([-0.05 0.5])
xlim([min(fICs)-0.05, max(fICs)+0.05])
xlabel('f_{IC}')
ylabel('Bias')
switch fittingtype
    case 'normal'
        title(['Rice (\sigma_0 = ' num2str(sigma0s(1)) ')' ])
    case 'adaptive'
        title(['Ratio (\sigma_0 = ' num2str(sigma0s(1)) ')' ])
end
grid on
ax = gca();
ax.FontSize=12;
legend;
f1.Position = [400, 200, 350, 400];

% STANDARD DEVIATION FIGURE
f2=figure;
for T2indx = 1:NT2
    T2 = T2s(T2indx);
    plot(fICs, sqrt(Variances(T2indx,:)), '-*', DisplayName = ['T2 = ' num2str(T2) ' ms'])
    hold on
end
ylim([-0.015 0.15])
xlim([min(fICs)-0.05, max(fICs)+0.05])
xlabel('f_{IC}')
ylabel('Standard deviation')
switch fittingtype
    case 'normal'
        title(['Rice (\sigma_0 = ' num2str(sigma0s(1)) ')' ])
    case 'adaptive'
        title(['Ratio (\sigma_0 = ' num2str(sigma0s(1)) ')' ])
end
grid on
ax = gca();
ax.FontSize=12;
% legend;
f2.Position = [400, 200, 350, 400];

