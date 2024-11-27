% Script to create new training data for MLP model training

% root folder
rootfolder = pwd;

% Number of training data samples
Nvoxel = 200;

% Noise 
noisetype = 'Rice';
sigma0s = [ 0.07];
T2s = [ 100];

% Training data folder
savedata = true;
TrainingDataFolder = fullfile(rootfolder, 'Simulation Experiment' , 'MLP', 'training data');


%% Protocols

% If multiple, define as a cell array of char arrays
modeltypes = {'Original VERDICT'};
schemenames = {'Original'};

% Saving sheme as mat file
schemesfolder = fullfile(rootfolder, 'Simulation Experiment','Schemes');
savescheme = true;


%% Cell sizes

% Cell Radius distribution R~Normal( muR, sigmaR)
muRmin = 6;
muRmax = 9;
sigmaRmin = 1;
sigmaRmax = 1;


%% Create training data

for sigma0 = sigma0s
    for T2 = T2s
        for indx = 1:length(modeltypes)
        
            modeltype = modeltypes{indx};
            schemename = schemenames{indx};
        
            disp([modeltype ' ' schemename])
        
            outputfolder = fullfile(TrainingDataFolder, modeltype, schemename, noisetype, ['T2_' num2str(T2)], ['sigma_' num2str(sigma0)] );
            
            createVERDICTdata( ...
                modeltype,...
                schemename,...
                Nvoxel = Nvoxel,...
                noisetype = noisetype,...
                sigma0 = sigma0,...
                T2 = T2,...
                randmuRs=[muRmin, muRmax],...
                randsigmaRs=[sigmaRmin, sigmaRmax],...
                schemesfolder=schemesfolder,...
                savedata=savedata,...
                outputfolder=outputfolder...
                );
        
        end
    end
end
