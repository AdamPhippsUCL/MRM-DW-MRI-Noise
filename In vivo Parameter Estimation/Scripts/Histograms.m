% Script to display histograms in figure 6

rootfolder = pwd;

% Patient
Patient = 'Patient_4';

% Load parameter maps
sigma0 = load(fullfile(rootfolder, 'In vivo Parameter Estimation', 'Outputs', Patient, 'sigma0.mat' )).sigma0;
T2 = load(fullfile(rootfolder, 'In vivo Parameter Estimation', 'Outputs', Patient, 'T2.mat' )).T2;


% Load ROI
ROI = load(fullfile(rootfolder , 'In vivo Parameter Estimation', 'Imaging Data' , 'ROIs', Patient, 'WP.mat')).img;
ROI = logical(ROI);

%% Create histograms

figfolder = fullfile(pwd, 'In vivo Parameter Estimation', 'Outputs', 'Figures', Patient);
mkdir(figfolder);

% sigma0
f1=figure;
sigma0bins = linspace(0, 0.08, 51);
h=histogram(sigma0(ROI), sigma0bins, FaceColor ="#0072BD");
hold on
f1.Position = [200,193,500,307];
ylim([0 1.05*max(h.Values)])
xlabel('\sigma_0')
ylabel('Counts')
ax = gca();
ax.FontSize = 12;
saveas(f1, fullfile(figfolder, 'sigma0hist.fig'));

% T2
f2=figure;
T2bins = linspace(0, 200, 51);
h=histogram(T2(ROI), T2bins, FaceColor=	"#D95319");
hold on
f2.Position = [200,200,500,300];
ylim([0 1.05*max(h.Values)])
xlabel('T2 (ms)')
ylabel('Counts')
ax = gca();
ax.FontSize = 12;
saveas(f2, fullfile(figfolder, 'T2hist.fig'));
