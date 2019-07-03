%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

%% Load in path data
load('proj.mat');

% temporary flag sets
proj.process.mri = 1;
proj.process.scr = 1;
proj.process.emg = 1;

proj.process.beta_mri_ex_id = 1;
proj.process.beta_scr_ex_id = 1;

%% Write out amended project structure
save('proj.mat','proj');
