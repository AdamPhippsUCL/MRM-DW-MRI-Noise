% Script to display example clibration curves for sigma0


% Root folder
rootfolder = pwd;

% Calibration folder
CalibrationCurveFolder = fullfile(rootfolder, 'In vivo Parameter Estimation', 'Functions', 'Sigma0 Calibration Curves');

% Define T2 values
T2s = [40, 50, 75, 100, 150, 200];


f=figure;
f.Position = [488 242 762.6000 420];

for T2 = T2s

    % Load calibration curve
    load(fullfile(CalibrationCurveFolder, ['T2 ' num2str(T2)], 'CalibrationCurve.mat'));

    % keys
    stds = keys(CalibrationCurve);

    % values
    sigma0s = values(CalibrationCurve);

    plot(stds, sigma0s, LineWidth = 1.5, DisplayName = ['T2 = ' num2str(T2)])
    hold on


end


ylim([0.005, 0.105])
xlim([0 0.3])
legend(Location="southeast")
xlabel('Normalized signal standard deviation')
ylabel('\sigma_0')
ax = gca();
ax.FontSize = 12;

% Save
saveas(f, fullfile(CalibrationCurveFolder, 'ExampleCurves.fig'))