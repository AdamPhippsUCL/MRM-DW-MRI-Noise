% Create MATLAB structure containing image details

% Define root folder
rootfolder = pwd;

% Folder to save structure
folder = fullfile(rootfolder, 'Phantom Experiments', 'Imaging Data');

% Initialise structure
ImageDetails = struct();


% Protocol names
ProtocolNames = {
    "b1500_Ex1",...
    "b2000_Ex1",...
    "b2100_Ex1",...
    "b2400_Ex1",...
    "b2500_Ex1",...
    "b2750_Ex1",...
    "b3000_Ex1",...
    "b3500_Ex1",...
    ...
    "b2500_Nb2_Ex2",...
    "b2500_Nb4_Ex2",...
    "b2500_Nb8_Ex2",...
    "b2500_Nb12_Ex2",...
    "b2500_Nb16_Ex2",...    
    };


% b values
bvalues = {
    1500,...
    2000,...
    2100,...
    2400,...
    2500,...
    2750,...
    3000,...
    3500,...
    ...
    2500,...
    2500,...
    2500,...
    2500,...
    2500,...
};


% TEs
TEs = {
    450,...
    200,...
    100,...
    234,...
    300,...
    450,...
    400,...
    150,...
    ...
    300,...
    300,...
    300,...
    300,...
    300,...
    };

% TRs
TRs = {
    4000,...
    4000,...
    4000,...
    4000,...
    4000,...
    4000,...
    4000,...
    4000,...
    ...
    4000,...
    4000,...
    4000,...
    4000,...
    4000,...
    };


% NSA
NSAs = {
    1,...
    1,...
    1,...
    1,...
    1,...
    1,...
    1,...
    1,...
    ...
    2,...
    2,...
    2,...
    2,...
    2,...
};


% Rav
Ravs = {
    3,...
    3,...
    3,...
    3,...
    3,...
    3,...
    3,...
    3,...
    ...
    3*1,...
    3*2,...
    3*4,...
    3*6,...
    3*8,...
    };


% Append all to structure
ImageDetails.ProtocolName = ProtocolNames;
ImageDetails.bvalue = bvalues;
ImageDetails.TE = TEs;
ImageDetails.TR = TRs;
ImageDetails.NSA = NSAs;
ImageDetails.Rav = Ravs;


% Save structure as mat file
save(fullfile(folder,'ImageDetails.mat'), 'ImageDetails')

