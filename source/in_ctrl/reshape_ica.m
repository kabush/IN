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

%% ----------------------------------------
%% ----------------------------------------
%% Configure the ACC mask

%% copy ica to tmp
eval(['! cp ',proj.path.atlas,'ray2013/ica70/maps/' ...
                    'thresh_zstatd70_17.nii.gz ',proj.path.code,'tmp/']);
eval(['! cp ',proj.path.atlas,'TT/TT_icbm452_orig.nii ',proj.path.code,'tmp/']);

%% rotate the ICA to match the rotation of our fMRI data
eval(['! 3dresample -orient RAI -prefix ',proj.path.code,'tmp/orient_thresh_zstatd70_17.nii.gz ' ...
      '-input ',proj.path.code,'tmp/thresh_zstatd70_17.nii.gz ']);

%% move to permanent storage
eval(['! mv ',proj.path.code,'tmp/orient_thresh_zstatd70_17.nii.gz ',proj.path.ctrl.in_ica]);
eval(['! mv ',proj.path.code,'tmp/TT_icbm452_orig.nii ',proj.path.ctrl.in_ica]);

%% clean-up
eval(['! rm ',proj.path.code,'tmp/*']);

%% ----------------------------------------
%% ----------------------------------------
%% Configure the State mask (for EVC)

%% copy ica to tmp
eval(['! cp ',proj.path.atlas,'ray2013/ica20/maps/*.gz ',proj.path.code,'tmp/']);

%% rotate the ICA to match the rotation of our fMRI data

% emotion/interoception
eval(['! 3dresample -orient RAI -prefix ',proj.path.code,'tmp/orient_thresh_zstatd20_1.nii.gz ' ...
      '-input ',proj.path.code,'tmp/thresh_zstat1.nii.gz ']);

eval(['! 3dresample -orient RAI -prefix ',proj.path.code,'tmp/orient_thresh_zstatd20_2.nii.gz ' ...
      '-input ',proj.path.code,'tmp/thresh_zstat2.nii.gz ']);

eval(['! 3dresample -orient RAI -prefix ',proj.path.code,'tmp/orient_thresh_zstatd20_3.nii.gz ' ...
      '-input ',proj.path.code,'tmp/thresh_zstat3.nii.gz ']);

eval(['! 3dresample -orient RAI -prefix ',proj.path.code,'tmp/orient_thresh_zstatd20_4.nii.gz ' ...
      '-input ',proj.path.code,'tmp/thresh_zstat4.nii.gz ']);

eval(['! 3dresample -orient RAI -prefix ',proj.path.code,'tmp/orient_thresh_zstatd20_5.nii.gz ' ...
      '-input ',proj.path.code,'tmp/thresh_zstat5.nii.gz ']);

% higher cognition
% ICAs: 13,14,15,16,17,18

%% move to permanent storage
eval(['! mv ',proj.path.code,'tmp/orient_thresh_zstatd20_1.nii.gz ',proj.path.ctrl.in_ica]);
eval(['! mv ',proj.path.code,'tmp/orient_thresh_zstatd20_2.nii.gz ',proj.path.ctrl.in_ica]);
eval(['! mv ',proj.path.code,'tmp/orient_thresh_zstatd20_3.nii.gz ',proj.path.ctrl.in_ica]);
eval(['! mv ',proj.path.code,'tmp/orient_thresh_zstatd20_4.nii.gz ',proj.path.ctrl.in_ica]);
eval(['! mv ',proj.path.code,'tmp/orient_thresh_zstatd20_5.nii.gz ',proj.path.ctrl.in_ica]);

%% clean-up
eval(['! rm ',proj.path.code,'tmp/*']);
