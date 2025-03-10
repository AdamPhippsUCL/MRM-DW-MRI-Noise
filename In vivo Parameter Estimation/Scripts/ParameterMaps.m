% Script to display parameter maps in figure 5

rootfolder = pwd;

% Patient
Patient = 'Patient_3';

% Load parameter maps
b0 = load(fullfile(rootfolder , 'In vivo Parameter Estimation', 'Imaging Data' , 'MAT', Patient, 'avgTE1.mat')).avgTE1;
sigma0 = load(fullfile(rootfolder, 'In vivo Parameter Estimation', 'Outputs', Patient, 'sigma0.mat' )).sigma0;
T2 = load(fullfile(rootfolder, 'In vivo Parameter Estimation', 'Outputs', Patient, 'T2.mat' )).T2;


% Load ROI
ROI = load(fullfile(rootfolder , 'In vivo Parameter Estimation', 'Imaging Data' , 'ROIs', Patient, 'WP.mat')).img;
ROI = double(ROI);

% Define slice and image region to display
switch Patient
    case 'Patient_1'
        slice = 12;
        dispinds = [27 82 23 78];
    case 'Patient_2'
        slice = 11;
        dispinds = [28 88 29 84];
    case 'Patient_3'
        slice = 8;
        dispinds = [27 79 30 82];
    case 'Patient_4'
        slice = 8;
        dispinds = [27 79 32 84];
end


%% Display images

figfolder = fullfile(pwd, 'In vivo Parameter Estimation', 'Outputs', 'Figures', Patient);
mkdir(figfolder);

f1=figure;
f1.Position = [300   200   409   300];
tiledlayout(1,1, "TileSpacing","compact");
nexttile;
b0ROI = b0.*ROI;
imshow(b0ROI(dispinds(3):dispinds(4),dispinds(1):dispinds(2), slice), []);
c=colorbar;
c.Label.String = 'b=0 signal';
ax = gca();
ax.FontSize = 12;
saveas(f1, fullfile(figfolder, 'b0.fig'));


f2=figure;
f2.Position = [300   200   409   300];
tiledlayout(1,1, "TileSpacing","compact");
nexttile;
sigma0ROI = sigma0.*ROI;
imshow(sigma0ROI(dispinds(3):dispinds(4),dispinds(1):dispinds(2), slice), [0.0 0.08]);
c=colorbar;
c.Label.String = '\sigma_0';
ax = gca();
ax.FontSize = 12;
saveas(f2, fullfile(figfolder, 'sigma0.fig'));

f3=figure;
f3.Position = [300   200   409   300];
tiledlayout(1,1, "TileSpacing","compact");
nexttile;
T2ROI = T2.*ROI;
imshow(T2ROI(dispinds(3):dispinds(4),dispinds(1):dispinds(2), slice), [0 200]);
c=colorbar;
c.Label.String = 'T2 (ms)';
ax = gca();
ax.FontSize = 12;
saveas(f3, fullfile(figfolder, 'T2.fig'));