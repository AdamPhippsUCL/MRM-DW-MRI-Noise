% Script to perform phantom experiment 2

% Define root folder
rootfolder = pwd;

% Imaging data folder
ImageDatafolder = fullfile(rootfolder, 'Phantom Experiments', 'Imaging Data');


%% Experiment Settings

% Noise model in fitting
noisetype = 'Ratio';

% Use denoised images?
usedenoised = true;

% Figure visibility
figvis = 'off';

% Define image names to use
ImageNames = {
    'b2500_Nb2_Ex2',...
    'b2500_Nb4_Ex2',...
    'b2500_Nb8_Ex2',...
    'b2500_Nb12_Ex2',...
    'b2500_Nb16_Ex2',...
    };

Nimg = length(ImageNames);

% Define ROI names to use
ROINames = {
    "ROI1",...
    "ROI2",...
    % "ROI3",...
    % "ROI4",...
    % "ROI5",...
    % "ROI6",...
    % "ROI7",...
    % "ROI8"
    };

NROI = length(ROINames);


%% Define folders

switch usedenoised
    case true
        imgtype = 'MAT DN';
    case false
        imgtype = 'MAT';
end

imagefolder = fullfile(ImageDatafolder, imgtype);
ROIfolder = fullfile(ImageDatafolder , 'ROIs');

% Load image and ROI details
load([ImageDatafolder '/ImageDetails.mat'])
load([ImageDatafolder '/ROIDetails.mat'])



%% Experiment

% Initial structure for results
FittingResults = struct();

for ROIindx = 1:NROI

    % Get ROI details
    ROIName = ROINames{ROIindx};
    whereROI = [ROIDetails.ROIName{:}] == ROIName;
    T2 = ROIDetails.T2{whereROI};
    D = ROIDetails.D{whereROI};


    for imgindx = 1:Nimg

        % Get image details
        ImageName = ImageNames{imgindx};
        whereimg = [ImageDetails.ProtocolName{:}] == ImageName;
        b = ImageDetails.bvalue{whereimg};
        TE = ImageDetails.TE{whereimg};
        TR = ImageDetails.TR{whereimg};
        NSA = ImageDetails.NSA{whereimg};
        Rav = ImageDetails.Rav{whereimg};

        disp(['ROI: ' char(ROIName) ', Image: ' char(ImageName)])
        

        % Load ROI
        ROIstruct = load(fullfile(ROIfolder , char(ImageName) , [char(ROIName) '.mat']));
        for name = fieldnames(ROIstruct)
            ROI = getfield(ROIstruct, name{1});
            ROI = (ROI==1);
        end

        % Load image
        img = double(load(fullfile(imagefolder , char(ImageName) , 'ImageArray.mat')).ImageArray);

        % Normalised image
        normimg = img(:,:,:,2)./img(:,:,:,1);

        % Remove nans and infinities
        normimg(isnan(normimg)) = 0;
        normimg(isinf(normimg)) = 0;

        % Extract ROI values
        ROIvalues = normimg(ROI);

        % Rescale if using DN images (rescale map from bias removal method)
        if strcmp(imgtype, 'MAT DN')
            Rescale = double(load(fullfile(imagefolder , char(ImageName) , 'Rescale.mat')).Rescale);
            ROIrescale = Rescale(ROI);
            shift = mean(ROIvalues)*(1-mean(ROIrescale));
            ROIvalues = ROIvalues - shift;
        end


        % Define histogram bins
        binmin = 0;
        binmax = 2;
        if max(ROIvalues)<1
            binmax = 1;
        end
        
        % Temporary histogram to get optimal bin spacing 
        ftemp=figure(Visible="off");
        h = histogram(ROIvalues, 50);
        binedges = h.BinEdges;
        binspacing = binedges(2)-binedges(1);
        close(ftemp)

        % Max bins 1000, min bins 100
        nbin =  min( [round((binmax-binmin)/binspacing), 1000] );
        if nbin<100
            nbin=100;
        end
        

        binedges = linspace(binmin, binmax, nbin+1);
        bincentres = (binedges(:,1:end-1) + binedges(:,2:end)) / 2 ;
        binspacing = bincentres(2)-bincentres(1);

        f=figure('visible',figvis);
        f.Position = [400, 200, 600, 400];
        h = histogram(ROIvalues, binedges, HandleVisibility='off');
        counts = h.Values;
        xlim([0 binmax])
        ylim([0 1.05*max(counts)]);


        % == Fitting noise model to histogram

        % First guesses
        sigma0guess = sqrt(2)*std(ROIvalues)/exp(TE/T2);
        fdguess = median(ROIvalues);

        beta0guess = [sigma0guess, fdguess];

        % Bounds
        Derr = 2e-4;
        sigma0min = 0.001;
        sigma0max = 0.4;
        lb = [sigma0min, exp(-(D+Derr)*b)];
        ub = [sigma0max, exp(-(D-Derr)*b)];

        % Apply fitting
        [coeffs, resnorm] = fitDistToHist( ...
            counts, ...
            bincentres, ...
            disttype=noisetype,...
            T2 = T2,...
            TE=TE, ...
            N0 = NSA, ...
            Nb = NSA*Rav,...
            beta0guess=beta0guess,...
            lb = lb, ...
            ub = ub);

        
        % Fitted parameters
        sigma0fit = coeffs(1);
        fdfit = coeffs(2);

        Dfit = (-1/b)*log(fdfit);


        % == Plot fitted distribution
        b0signal = exp(-TE/T2);
        bsignal = fdfit*b0signal;

        switch noisetype
            case 'Ratio'
                [dist, signals] = RatioDist(b0signal, bsignal, sigma0fit, N0 = NSA, Nb = NSA*Rav, zs = bincentres);
            case 'Rice'
                [dist, signals] = RiceDist(b0signal, bsignal, sigma0fit, zs = bincentres);
        end
        
        hold on

        xlabel('Normalized signal')
        ylabel('Counts')
        switch noisetype
            case 'Ratio'
                plot(signals, dist*binspacing*sum(counts), LineWidth = 2, DisplayName = [noisetype ' distribution'], color = "#D95319");
            case 'Rice'
                plot(signals, dist*binspacing*sum(counts), LineWidth = 2, DisplayName = [noisetype ' distribution'], color = "#7E2F8E");
        end

        ax = gca();
        ax.FontSize = 12;
        legend;
        if strcmp(figvis, 'on')
            pause(1)
        end
        close(f);

        % Fill in results
        FittingResults((ROIindx-1)*Nimg + imgindx).ImageName = ImageName;
        FittingResults((ROIindx-1)*Nimg + imgindx).ROIName = ROIName;        
        FittingResults((ROIindx-1)*Nimg + imgindx).sigma0fit = sigma0fit;
        FittingResults((ROIindx-1)*Nimg + imgindx).T2 = T2;
        FittingResults((ROIindx-1)*Nimg + imgindx).D = D;        
        FittingResults((ROIindx-1)*Nimg + imgindx).Dfit = Dfit;
        FittingResults((ROIindx-1)*Nimg + imgindx).resnorm = resnorm;

    end
end


%% Plot results

colordict = dictionary([1,2,3,4,5,6,7,8], [	"#0072BD", 	"#D95319",	"#EDB120", 	"#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F", "#e981cd"]);
concentrations = dictionary([1,2,3,4,5,6,7,8], ["50% PVP", "40% PVP", "30% PVP", "20% PVP", "10% PVP", "5% PVP", "2.5% PVP", "1mM NiCl_2"]);

Nbs = zeros(1,Nimg);
for imgindx = 1:Nimg
    % Get image details
    ImageName = ImageNames{imgindx};
    whereimg = [ImageDetails.ProtocolName{:}] == ImageName;
    NSA = ImageDetails.NSA{whereimg};
    Rav = ImageDetails.Rav{whereimg};
    Nbs(imgindx)=NSA*Rav;
end


fig1 = figure;
fig1.Position = [300   200   400   400];

sigma0s = zeros(NROI, Nimg);

% Scatter plot for each ROI
for ROIindx = 1:NROI

    % Get ROI details
    ROIName = ROINames{ROIindx};
    whereROI = [ROIDetails.ROIName{:}] == ROIName;
    T2 = ROIDetails.T2{whereROI};
    D = ROIDetails.D{whereROI};

    % ROI bools
    ROIbools = ([FittingResults.ROIName] == ROIName);

    % Extract sigma0 measurements
    sigma0s(ROIindx, :) = [FittingResults(ROIbools).sigma0fit];


    plot(1./sqrt(Nbs), sigma0s(ROIindx,:)/mean(sigma0s(ROIindx, :)), '*-', DisplayName = concentrations(ROIindx), Color = colordict(ROIindx))
    hold on
    xlabel('N_{b>0}')
    ylabel('\sigma_0 / mean(\sigma_0)')
    if ~usedenoised
        title(noisetype)
    else
        title([char(noisetype) ' (DN)'])
    end
    legend(NumColumns = 2)
    ax = gca();
    ax.FontSize = 12;

end

%% Save figures and fitting results

outputfolder = fullfile(rootfolder, 'Phantom Experiments', 'Outputs');
dt = char(datetime());
dt = strrep(dt, ':', '-');
outf = fullfile(outputfolder , dt);
figfolder = fullfile(outf ,'figures');
mkdir(figfolder);


% Meta information
Meta = struct();
Meta.imagefolder = imgtype;
Meta.ImageNames = ImageNames;
Meta.ROINames = ROINames;
Meta.NoiseType = noisetype;
Meta.Derr = Derr;
Meta.sigma0range = [sigma0min, sigma0max];

save([outf '/Meta.mat'], "Meta");
save([outf '/FittingResults.mat'], "FittingResults");

% Save figures
saveas(fig1, [char(figfolder) '/sigma0.fig'])
clear;