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

%% find ica values greater than threshold to achieve largest single
%% ROI (>5.6734)
eval(['! 3dcalc -a ',proj.path.code,'tmp/orient_thresh_zstatd70_17.nii.gz ' ...
                    '-expr  ''ispositive(a-5.6734)'' -prefix ' ...
                    ,proj.path.code,'tmp/frc_orient_thresh_zstatd70_17.nii.gz']);

%% change from 1x1x1 to 3x3x3 voxel sizes
eval(['! 3dfractionize -template ',proj.path.mri.gm_mask,'group_gm_mask.nii ' ...
      '-input ',proj.path.code,'tmp/frc_orient_thresh_zstatd70_17.nii.gz ' ...
      '-prefix ',proj.path.code,'tmp/' ...
      'frc_orient_thresh_zstatd70_17_3x3x3.nii.gz -clip .2']);

%% create mask 
eval(['! 3dcalc -a ',proj.path.code,'tmp/frc_orient_thresh_zstatd70_17_3x3x3.nii.gz ' ...
      '-expr  ''bool(a)'' -prefix ' ...
      ,proj.path.code,'tmp/sng_orient_thresh_zstatd70_17_3x3x3.nii.gz']);


%% move to permanent storage
eval(['! mv ',proj.path.code,'tmp/sng_orient_thresh_zstatd70_17_3x3x3.nii.gz ',proj.path.ctrl.in_ica]);
eval(['! mv ',proj.path.code,'tmp/TT_icbm452_orig.nii ',proj.path.ctrl.in_ica]);

%% clean-up
eval(['! rm ',proj.path.code,'tmp/*']);

%% ----------------------------------------
%% ----------------------------------------
%% Configure the State mask (for EVC)

%% copy ica to tmp
eval(['! cp ',proj.path.atlas,'ray2013/ica20/maps/*.gz ',proj.path.code,'tmp/']);

%% rotate the ICA to match the rotation of our fMRI data

%% ----------------------------------------
%% EMOTION/INTEROCEPTION REGIONS
%%
for i=1:5

    %% orient ICA to RAI like the rest of the data
    eval(['! 3dresample -orient RAI -prefix ',proj.path.code,'tmp/' ...
                        'orient_thresh_zstatd20_',num2str(i),'.nii.gz ' ...
                        '-input ',proj.path.code,'tmp/thresh_zstat',num2str(i),'.nii.gz ']);
    
    %% change from 1x1x1 to 3x3x3 voxel sizes
    eval(['! 3dfractionize -template ',proj.path.mri.gm_mask,'group_gm_mask.nii ' ...
          '-input ',proj.path.code,'tmp/orient_thresh_zstatd20_',num2str(i),'.nii.gz ' ...
          '-prefix ',proj.path.code,'tmp/' ...
          'orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz -clip .2']);

    %% create mask
    eval(['! 3dcalc -a ',proj.path.code,'tmp/orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz ' ...
          '-expr  ''bool(a)'' -prefix ' ...
          ,proj.path.code,'tmp/sng_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz']);

end

%% ----------------------------------------
%% HIGHER COGNITION REGIONS
%%
% TBD ICAs: 13,14,15,16,17,18

%% move to permanent storage
eval(['! mv ',proj.path.code,'tmp/sng_orient_thresh_zstatd20_1_3x3x3.nii.gz ',proj.path.ctrl.in_ica]);
eval(['! mv ',proj.path.code,'tmp/sng_orient_thresh_zstatd20_2_3x3x3.nii.gz ',proj.path.ctrl.in_ica]);
eval(['! mv ',proj.path.code,'tmp/sng_orient_thresh_zstatd20_3_3x3x3.nii.gz ',proj.path.ctrl.in_ica]);
eval(['! mv ',proj.path.code,'tmp/sng_orient_thresh_zstatd20_4_3x3x3.nii.gz ',proj.path.ctrl.in_ica]);
eval(['! mv ',proj.path.code,'tmp/sng_orient_thresh_zstatd20_5_3x3x3.nii.gz ',proj.path.ctrl.in_ica]);

%% clean-up
eval(['! rm ',proj.path.code,'tmp/*']);
