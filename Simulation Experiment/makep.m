
% RUN SCRIPT TO TURN ALL FUNCTIONS IN 'FUNCTIONS_M' FOLDER INTO PCODE
% FUNCTIONS IN NEW 'FUNCTIONS' FOLDER.

rootfolder = pwd;

thisfolder = fullfile(rootfolder, 'Simulation Experiment');

% Copy Functions folder
copyfile(fullfile(thisfolder, 'Functions_M'), fullfile(thisfolder, 'Functions_P') );

% Find all function .m files in 'Functions' folder
FileList = dir(fullfile(thisfolder, 'Functions_P', '**', '*.m'));

for findx = 1:length(FileList)

    disp(FileList(findx).name);

    fname = fullfile(FileList(findx).folder, FileList(findx).name);

    pcode(fname, "-inplace");

    delete(fname)


end