% Initializes the L-TUTAM
% Lagrangian Tampere University of Technology Aerosol Model
% Miska Olin 4.5.2016

disp('L-TUTAM starting...')
if ~exist('p','var')
    disp('Loading interpolation tables for the PL distribution ...')
    load taulukot
end

% Load the default setup
defaultSetupFile;

% Asking for own setup file
disp('Select the setup file')
filename = 'defaultSetupFile.m';
pathname = '';
[filename,pathname] = uigetfile('*.m','Select the setup file',strcat(pathname,filename));


% Runnig your own setup file
run(strcat(pathname,filename));

clear filename pathname

disp(p)




