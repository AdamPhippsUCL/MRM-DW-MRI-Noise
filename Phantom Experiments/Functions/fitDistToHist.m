function [coeffs, resnorm] = fitDistToHist(counts, bincenters, params, opts)

arguments

    counts % Array of counts
    bincenters % Array of bin centers

    % Distribution parameters 
    params.sigma0 = NaN
    params.fd = NaN
    params.T2 
    params.TE = 50 % Specify at [TE1, TE2] for DTIRatio
    params.N0 = 1
    params.Nb = 1

    opts.disttype 
    opts.beta0guess 

    opts.lb 
    opts.ub 

end


%% Normalise histogram counts

counts = counts/sum(counts);


%% Fit

% Assign variable to base (for access in functions)
assignin('base', 'current_params', params)

options = optimoptions('lsqcurvefit');


switch opts.disttype

    case 'Ratio'
        [coeffs, resnorm] = lsqcurvefit(@evalRatioDist, opts.beta0guess, bincenters, counts, opts.lb, opts.ub, options);

    case 'Rice'
        [coeffs, resnorm] = lsqcurvefit(@evalRiceDist, opts.beta0guess, bincenters, counts, opts.lb, opts.ub , options);  

end


end


% Function to fit Ratio distribution
function binfreqs = evalRatioDist(b, x)

arguments
    b % [sigma0, fd] value
    x % bin centers
end

% Get parameters
params = evalin('base', 'current_params');
T2 = params.T2;
TE = params.TE;
N0 = params.N0;
Nb = params.Nb;
sigma0 = b(1);
fd = b(2);

% b0 signal
b0signal = exp(-TE/T2);

% b signal
bsignal = fd*b0signal;

% Bin center
bincentres = x;

% Bin spacings
binmin = min(x);
binmax = max(x);
binspacing = x(2)-x(1);

% Generate pdf over bin centers
pdfvals = RatioDist( ...
    b0signal, ...
    bsignal, ...
    sigma0, ...
    N0=N0, ...
    Nb=Nb, ...
    zs = bincentres,...
    ys = bincentres,...
    ymin = binmin,...
    ymax = binmax,...
    dy = binspacing);

% Evaluate bin frequencies
binfreqs = pdfvals*binspacing;

binfreqs = reshape(binfreqs, size(x));

end


% % Function to fit Rice distribution
function binfreqs = evalRiceDist(b, x)

arguments
    b % [sigma0, fd] values
    x % bin centers
end

% Get parameters
params = evalin('base', 'current_params');
sigma0 = b(1);
fd = b(2);

% Bin center
bincentres = x;

% Bin spacings
binmin = min(x);
binmax = max(x);
binspacing = x(2)-x(1);

% Generate pdf over bin centers
pdfvals = RiceDist(1, fd, sigma0, zs = bincentres);

% Evaluate bin frequencies
binfreqs = pdfvals*binspacing;

binfreqs = reshape(binfreqs, size(x));

end

