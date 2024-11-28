% Script to generate calibration curves for sigma0

% Root folder
rootfolder = pwd;


%% Simulation

% Define echo times
TE1 = 50;
TE2 = 75;

% Define T2 values
T2s = linspace(40, 240, 41);

% Define sigma0 values
sigma0s = linspace(0.01, 0.1, 91);

% Number of signal samples
Nsample = 10000;

% Initialise array
simulatedstds = zeros(size(sigma0s));

for T2 = T2s
    for sigma0indx = 1:length(sigma0s)
    
        sigma0 = sigma0s(sigma0indx);
        
        % Generate ratio distribution
        A1 = exp(-TE1/T2);
        A2 = exp(-TE2/T2);
        [dist, signals] = RatioDist(A1, A2, sigma0); 

        % Take samples
        samples = zeros(1, Nsample);
        for sampleindx = 1:Nsample
            samples(sampleindx)=sampleDistribution(dist, signals);
        end

        % Append standard deviation
        simulatedstds(sigma0indx) = std(samples);

    end

    % figure;
    % plot(simulatedstds, sigma0s)
    % hold on

    % Savitzky-Golay filter
    filtwidth = 21;
    filtorder = 2;
    simulatedstds = sgolayfilt(simulatedstds, filtorder, filtwidth);

    % plot(simulatedstds, sigma0s)


    % Create dictionary
    CalibrationCurve = dictionary(simulatedstds, sigma0s);


    % Meta information
    Meta = struct();
    Meta.Nsample = Nsample;
    Meta.sigma0s = sigma0s;
    Meta.filtwidth = filtwidth;
    Meta.filtorder = filtorder;

    % Save dictionary
    folder = fullfile(rootfolder, 'In vivo Parameter Estimation', 'Functions', 'Sigma0 Calibration Curves', ['T2 ' num2str(T2)]);
    mkdir(folder);
    save(fullfile(folder, 'CalibrationCurve.mat'), "CalibrationCurve");
    save(fullfile(folder, 'Meta.mat'), "Meta");

end


% Save T2 values
save(fullfile(rootfolder, 'In vivo Parameter Estimation', 'Functions', 'Sigma0 Calibration Curves', 'T2s.mat'), "T2s");

