% MATLAB function to train an MLP to fit VERDICT model
% Function calls python script to train and export MLP model
function trainMLP(datafolder, modelfolder, opts)


arguments
    
    datafolder % folder containing training data
    modelfolder % folder to save MLP models

    % MLP Architecture
    opts.Nlayer = 3
    opts.Nnode = 150
    
    % MLP training
    opts.batchsize = 100

    % Epochs
    opts.maxepochs = 150


end


% == Load training data

load([datafolder '/Meta.mat']);
load([datafolder '/params.mat']);
load([datafolder '/signals.mat']);

% Scale data
% [scaledSignals, minSignal, maxSignal] = minmaxscale(Signals);
% signalscaling = [minSignal, maxSignal];
scaledSignals=Signals;
signalscaling = [0, 1];

% [scaledParams, minParams, maxParams] = minmaxscale(params);
% paramscaling = [minParams, maxParams];
scaledParams=params;
paramscaling  = [0, 1];


% == Configure network

% Sizes of layers
Ninput = size(Signals, 2);
Noutput = size(params, 2);

% Define layers
HiddenLayer = [fullyConnectedLayer(opts.Nnode), reluLayer];
HiddenLayers = repmat(HiddenLayer, 1, opts.Nlayer);

layers = [ ...
    sequenceInputLayer(Ninput),...
    HiddenLayers,...
    fullyConnectedLayer(Noutput)
    ] ;

net = dlnetwork(layers) ;

% == Training options
options = trainingOptions( ...
    "adam", ...
    MaxEpochs=opts.maxepochs, ...
    Verbose=false, ...
    Plots="training-progress", ...
    Metrics="rmse",...
    LearnRateSchedule='none');


lossFcn = "l2loss" ; % Loss between true and predicted LWF 

% == Apply training
[mlp, info] = trainnet(scaledSignals, scaledParams, net, lossFcn, options);

% Close figure
delete(findall(0));


% == Save MLP and scaling information

% Create output folder
outputfolder = modelfolder;
mkdir(outputfolder);

save([outputfolder '/mlp.mat'], "mlp");
save([outputfolder '/signalscaling.mat'], "signalscaling");
save([outputfolder '/paramscaling.mat'], "paramscaling");
     




%% META data

FitMeta = struct();
FitMeta.train_complete_time = datetime();
FitMeta.TrainingDataMeta = Meta;
FitMeta.Nlayer = opts.Nlayer;
FitMeta.Nnode = opts.Nnode;
FitMeta.batchsize = opts.batchsize;

save([modelfolder '/Meta.mat'], 'FitMeta');




end


% Min Max scaler function
function [scaledData, minVal, maxVal] = minmaxscale(data)
    % Get the minimum and maximum values
    minVal = min(data(:));
    maxVal = max(data(:));
    
    % Scale the data to the range [0, 1]
    scaledData = (data - minVal) ./ (maxVal - minVal);
end