% Script to train MLP models

% root folder
rootfolder = pwd;

% Noise 
noisetype = 'Rice';
sigma0s = [ 0.02, 0.04, 0.06];
T2s = [ 10000];


% Specify protocols 
modeltypes = {'Original VERDICT'} ;
schemenames = {'Original'}; 

TrainingDataFolder = fullfile(rootfolder, 'Simulation Experiment', 'MLP', 'training data');
ModelFolder = fullfile(rootfolder, 'Simulation Experiment', 'MLP', 'models');

for sigma0 = sigma0s
    for T2 = T2s

        for indx = 1:length(modeltypes)
            
            modeltype = modeltypes{indx};
            schemename = schemenames{indx};

            datafolder = fullfile(TrainingDataFolder , modeltype , schemename , noisetype, ['T2_' num2str(T2)], ['sigma_' num2str(sigma0)]);
            modelfolder = fullfile(ModelFolder, modeltype, schemename, noisetype, ['T2_' num2str(T2)],  ['sigma_' num2str(sigma0)] );
        
            trainMLP( ...
                datafolder,...
                modelfolder...
                );
            
        end

    end
end


for indx = 1:length(modeltypes)
    modeltype = modeltypes{indx};
    schemename = schemenames{indx};
    save(fullfile(ModelFolder, modeltype, schemename, noisetype, 'T2s.mat' ), "T2s");
    save(fullfile(ModelFolder, modeltype, schemename, noisetype, 'sigma0s.mat' ), "sigma0s");
end