% Function to call MATLAB MLP networks

function x = callMLP(signals, modelfolder)

arguments
    signals % Input signals
    modelfolder 
end


% Load MLP, scaling information, and meta data
modelfolder = char(modelfolder);
load([modelfolder '/Meta.mat']);
load([modelfolder '/mlp.mat']);
load([modelfolder '/paramscaling.mat']);
load([modelfolder '/signalscaling.mat']);


% Apply scaling to signals
scaledsignals = apply_scaling(signals, signalscaling);

% Call MLP
pred = minibatchpredict(mlp, scaledsignals);

% Apply inverse scaling
x = apply_inv_scaling(pred, paramscaling);


end


% Min Max scaler function
function scaledData = apply_scaling(data, scaling)

    N = length(scaling);
    minVal = scaling(1:N/2);
    maxVal = scaling(N/2+1:end);
    % Scale the data to the range [0, 1]
    scaledData = (data - minVal) ./ (maxVal - minVal);
end


% Inverse Min Max scaler function
function data = apply_inv_scaling(Data, scaling)
    
    N = length(scaling);
    minVal = scaling(1:N/2);
    maxVal = scaling(N/2+1:end);  

    data = Data*(maxVal-minVal) + minVal;
end

