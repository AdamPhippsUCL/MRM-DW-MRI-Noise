% Script to define and make scheme

% Define root folder
rootfolder = pwd;


%% Define sequence parameters

% Scheme: [delta,
%          DELTA, 
%          b, 
%          TE, 
%          TR, 
%          NSA,  (Number of signal averages for b=0 image)
%          Rav   (Additional signal averaging factor for b>0 image)
%           ]  

schemename = 'Original';

V01 = [1, 2, 0, 87, 4000, 6, 1];
V1 = [20, 43.4, 3000, 87, 4000, 6, 3*3]; % Rav = (3 gradient directions)*(3 Avg high b factor)
V02 = [1, 2, 0, 75, 4341, 6, 1];
V2 = [14, 37.4, 2000, 75, 4000, 6, 3*3]; % Rav = (3 gradient directions)*(3 Avg high b factor)
V03 = [1,2,0,94,2000,6,1];
V3 = [23.5, 46.9, 1500, 94, 4000, 6, 3*3]; % Rav = (3 gradient directions)*(3 Avg high b factor)
V04 = [1,2,0,68,2000,6,1];
V4 = [10.5, 33.9,500 68, 4000, 6, 3*2]; % Rav = (3 gradient directions)*(2 Avg high b factor)
V05 = [1,2,0,54,2000,6, 1];
V5 = [3.5, 26.9,90, 54, 4000, 6, 3*1]; % Rav = (3 gradient directions)*(1 Avg high b factor)

Vs = [...
    V01; V1;...
    V02; V2;...
    V03; V3;...
    V04; V4;...
    V05; V5;...            
    ];


%% Build scheme
scheme = BuildScheme(Vs, schemename);


%% Save scheme
schemefolder = fullfile(rootfolder, 'Simulation Experiment', 'Schemes');
save( fullfile(schemefolder, [schemename '.mat'] ) , 'scheme')
