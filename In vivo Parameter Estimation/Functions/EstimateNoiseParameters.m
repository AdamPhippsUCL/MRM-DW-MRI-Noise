% MATLAB function to perform in vivo estimation of Ratio noise model parameters

function [sigma0, T2] = EstimateNoiseParameters(IMG1, IMG2, TEvec, opts)

    arguments
        IMG1 % Image array for TE=TE1 [nx, ny, nz, nrep]
        IMG2 % Image array for TE=TE2 [nx, ny, nz, nrep]
        TEvec % vector of TE values [TE1, TE2]  

        opts.patchsize = [1 1] % If estimation done with patches
        opts.maskprcnt = 75
    end



% Initial check
if ~all(size(IMG1) == size(IMG2))
    error('Images have different sizes')
end

% Get image size
[Nx, Ny, Nz, Nrep] = size(IMG1);

% Echo times
TE1 = TEvec(1);
TE2 = TEvec(2);

%% Create stack of normalized images

RatioImageStack = zeros(Nx, Ny, Nz, Nrep^2);

for repindx1 = 1:Nrep

    img1 = IMG1(:,:,:,repindx1);

    for repindx2 = 1:Nrep

        img2 = IMG2(:,:,:,repindx2);

        % Remove zeros
        img1(img1==0) = eps;
        img2(img2==0) = eps;
    
        % Calculate ratio
        ratioimg = double(img2)./double(img1);
        
        RatioImageStack(:,:,:, (repindx1-1)*Nrep+repindx2 ) = ratioimg;

    end

end


%% Parameter estimation

T2 = zeros([Nx, Ny, Nz]);
sigma0 = zeros([Nx, Ny, Nz]);

% Patch width
pwx = ceil((opts.patchsize(1)-1)/2);
pwy = ceil((opts.patchsize(2)-1)/2);

% Masking
mskpct = opts.maskprcnt/100;

for zindx = 1:Nz
    for yindx = pwy + Ny*(1-mskpct)/2:Ny*(1+mskpct)/2-pwy
        for xindx = pwx + Nx*(1-mskpct)/2:Nx*(1+mskpct)/2-pwx

            disp(['Voxel Indices: ' num2str(xindx) ', ' num2str(yindx), ', ' num2str(zindx)])

            normvals = RatioImageStack(xindx-pwx:xindx+pwx, yindx-pwy:yindx+pwy, zindx, :);
            normvals  = normvals(:);

            % Define outlier range
            iqr = (prctile(normvals, 75) - prctile(normvals, 25));
            up = median(normvals) + 1.5*iqr;
            down = median(normvals) - 1.5*iqr;
            
            % Remove outliers
            outlierbools = or(normvals>up, normvals<down);
            normvals(outlierbools) = [];


            % == Fitting

            % T2
            T2guess = -(TE2-TE1)/log(median(normvals));
            T2guess(T2guess<25) = 25;
            T2guess(T2guess>250) = 250;

            % sigma0
            rTE = (2*(TE1^2))/(TE1+TE2);
            sigma0guess = std(normvals)/(sqrt(2)*exp(rTE/T2guess));

            sigma0fit = sigma0guess;
            T2fit = T2guess;

            % Append to arrays
            T2(xindx, yindx, zindx) = T2fit;
            sigma0(xindx, yindx, zindx) = sigma0fit;


        end
    end
end

end