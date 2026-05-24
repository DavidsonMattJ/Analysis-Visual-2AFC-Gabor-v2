% set_myOnlineDirectories_vGaborv2


%derive homedir from files location:
homedir = fileparts(mfilename('fullpath')); % check, not working?


%UTS mac 
cd('/Users/164376/Documents/GitHub/Analysis-Visual-2AFC-Gabor-v2/');
%%
homedir = pwd;
addpath(genpath([homedir filesep 'Analysis']));
addpath(genpath([homedir filesep 'Raw_Data']));
addpath(genpath([homedir filesep 'Processed_Data']));
addpath(genpath([homedir filesep 'Figures']));

datadir=  [homedir filesep 'Raw_Data'];
procdatadir= [homedir filesep 'Processed_Data'];
figdir = [homedir filesep 'Figures'];