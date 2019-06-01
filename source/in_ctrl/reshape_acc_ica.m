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

logger(['*******************************************'],proj.path.logfile);
logger([' Reshape Ray 2013 ACC ICA to match our data'],proj.path.logfile);
logger(['*******************************************'],proj.path.logfile);

%% ----------------------------------------
%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.ctrl.in_ica]);
    eval(['! rm -rf ',proj.path.ctrl.in_ica]);
    disp(['Creating ',proj.path.ctrl.in_ica]);
    eval(['! mkdir ',proj.path.ctrl.in_ica]);
end

%% copy ica to tmp
eval(['! cp ',proj.path.atlas,'ray2013/ica70/maps/' ...
                    'thresh_zstatd70_17.nii.gz ',proj.path.code,'tmp/']);
eval(['! cp ',proj.path.atlas,'TT/TT_icbm452_orig.nii ',proj.path.code,'tmp/']);

%% rotate the ICA to match the rotation of our fMRI data
eval(['! 3dresample -orient RAI -prefix ',proj.path.code,'tmp/orient_thresh_zstatd70_17.nii.gz ' ...
      '-input ',proj.path.code,'tmp/thresh_zstatd70_17.nii.gz ']);

%% move to permanent storaage
eval(['! mv ',proj.path.code,'tmp/orient_thresh_zstatd70_17.nii.gz ',proj.path.ctrl.in_ica]);
eval(['! mv ',proj.path.code,'tmp/TT_icbm452_orig.nii ',proj.path.ctrl.in_ica]);

%% clean-up
eval(['! rm ',proj.path.code,'tmp/*']);