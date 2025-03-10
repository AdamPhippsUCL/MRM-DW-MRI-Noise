% MATLAB script to run Ratio noise model parameter estimation

% Root folder
rootfolder = pwd;

% Patient
Patient = 'Patient_4';

% Load images
IMG1 = load(fullfile(rootfolder, 'In vivo Parameter Estimation', 'Imaging Data', 'MAT', Patient, 'TE1.mat')).TE1;
IMG2 = load(fullfile(rootfolder, 'In vivo Parameter Estimation', 'Imaging Data', 'MAT', Patient, 'TE2.mat')).TE2;

% Define TE vector
TEvec = [50, 125];


% Calibration curve folder 
CalibrationCurveFolder = fullfile(rootfolder, 'In vivo Parameter Estimation', 'Functions', 'Sigma0 Calibration Curves');

% Run parameter estimation
[sigma0, T2] = EstimateNoiseParameters(IMG1, IMG2, TEvec, CalibrationCurveFolder=CalibrationCurveFolder);


% Slice wise smoothing
sigmafilt = 1;
for zindx = 1:size(sigma0, 3)
    sigma0(:,:,zindx) = imgaussfilt(sigma0(:,:,zindx), sigmafilt);
    T2(:,:,zindx) = imgaussfilt(T2(:,:,zindx), sigmafilt);
end


%% Save parameter maps

mkdir(fullfile(rootfolder, 'In vivo Parameter Estimation', 'Outputs', Patient));
save(fullfile(rootfolder, 'In vivo Parameter Estimation', 'Outputs', Patient, 'sigma0.mat'), "sigma0");
save(fullfile(rootfolder, 'In vivo Parameter Estimation', 'Outputs', Patient, 'T2.mat'), "T2");