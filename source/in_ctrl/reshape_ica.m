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

% logger(['*******************************************'],proj.path.logfile);
% logger([' Reshape Ray 2013 ACC ICA to match our data'],proj.path.logfile);
% logger(['*******************************************'],proj.path.logfile);
% 
logger(['*******************************************'],proj.path.logfile);
logger([' Reshape masks and Ray 2013 ICAs to match our data'],proj.path.logfile);
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
%% Configure the mPFC (i.e., ACC) mask

%% copy ica to tmp
eval(['! cp ',proj.path.atlas,'Vega2016/mfc_mask.nii.gz ',proj.path.code,'tmp/']);
eval(['! cp ',proj.path.atlas,'TT/TT_icbm452_orig.nii ',proj.path.code,'tmp/']);

%% rotate the mask to match the rotation of our fMRI data
eval(['! 3dresample -orient RAI  -prefix ' ...
      '',proj.path.code,'tmp/orient_mfc_mask.nii.gz ' '-input ' ...
      '',proj.path.code,'tmp/mfc_mask.nii.gz']);

%% change from 2x2x2 to 3x3x3 voxel sizes
eval(['! 3dfractionize -template ',proj.path.mri.gm_mask,'group_gm_mask.nii ' ...
      '-input ',proj.path.code,'tmp/orient_mfc_mask.nii.gz ' ...
      '-prefix ',proj.path.code,'tmp/' ...
      'orient_mfc_mask_3x3x3.nii.gz -clip .2']);

%% create mask 
eval(['! 3dcalc -a ',proj.path.code,'tmp/orient_mfc_mask_3x3x3.nii.gz ' ...
      '-expr  ''bool(a)'' -prefix ' ...
      ,proj.path.code,'tmp/sng_orient_mfc_mask_3x3x3.nii.gz']);

%% ----------------------------------------

%% widen mask by 1st voxel
eval(['! 3dcalc -a ',proj.path.code,'tmp/sng_orient_mfc_mask_3x3x3.nii.gz ' ...
      '-b a+i -c a-i -d a+j -e a-j -f a+k -g a-k ' ...
      '-expr  ''bool(a+b+c+d+e+f+g)'' -prefix ' ...
      ,proj.path.code,'tmp/erode_sng_orient_mfc_mask_3x3x3.nii.gz']);

%% widen mask by 2nd voxel
eval(['! 3dcalc -a ',proj.path.code,'tmp/erode_sng_orient_mfc_mask_3x3x3.nii.gz ' ...
      '-b a+i -c a-i -d a+j -e a-j -f a+k -g a-k ' ...
      '-expr  ''bool(a+b+c+d+e+f+g)'' -prefix ' ...
      ,proj.path.code,'tmp/erode2_sng_orient_mfc_mask_3x3x3.nii.gz']);

%% ----------------------------------------
%% ----------------------------------------
%% Configure the State mask (for EVC)

%% copy ica to tmp
eval(['! cp ',proj.path.atlas,'ray2013/ica20/maps/*.gz ',proj.path.code,'tmp/']);

%% ----------------------------------------
%% As per Ray et al., (2013)
%% EMOTION/INTEROCEPTION REGIONS: ICs 1-5
%% SENSORY/MOTOR: ICs 6-12
%% HIGHER COGNITION: ICs 13-18
%%
ica_seq = proj.param.ctrl.ica_ids;

for i = 1:numel(ica_seq) % generate all 18 and use only first 5 (emotion)

    ica = ica_seq(i);

    %% orient ICA to RAI like the rest of the data
    eval(['! 3dresample -orient RAI -prefix ',proj.path.code,'tmp/' ...
          'orient_thresh_zstatd20_',num2str(ica),'.nii.gz ' ...
          '-input ',proj.path.code,'tmp/thresh_zstat',num2str(ica),'.nii.gz ']);
    
    %% change from 1x1x1 to 3x3x3 voxel sizes
    eval(['! 3dfractionize -template ',proj.path.mri.gm_mask,'group_gm_mask.nii ' ...
          '-input ',proj.path.code,'tmp/orient_thresh_zstatd20_',num2str(ica),'.nii.gz ' ...
          '-prefix ',proj.path.code,'tmp/' ...
          'orient_thresh_zstatd20_',num2str(ica),'_3x3x3.nii.gz -clip .2']);

    %% create full-ICA Mask
    eval(['! 3dcalc -a ',proj.path.code,'tmp/orient_thresh_zstatd20_',num2str(ica),'_3x3x3.nii.gz ' ...
          '-expr  ''bool(a)'' -prefix ' ...
          ,proj.path.code,'tmp/full_orient_thresh_zstatd20_',num2str(ica),'_3x3x3.nii.gz']);
    
    %% subtract out dACC from Mask
    eval(['! 3dcalc -a ',proj.path.code,'tmp/full_orient_thresh_zstatd20_',num2str(ica),'_3x3x3.nii.gz ' ...
          '-b ',proj.path.code,'tmp/erode2_sng_orient_mfc_mask_3x3x3.nii.gz ' ...
          '-expr  ''a-b'' -prefix ' ...
          ,proj.path.code,'tmp/diff_orient_thresh_zstatd20_',num2str(ica),'_3x3x3.nii.gz']);
    
    %% keep only the positive
    eval(['! 3dcalc -a ',proj.path.code,'tmp/diff_orient_thresh_zstatd20_',num2str(ica),'_3x3x3.nii.gz ' ...
          '-expr  ''ispositive(a)'' -prefix ' ...
          ,proj.path.code,'tmp/sng_orient_thresh_zstatd20_',num2str(ica),'_3x3x3.nii.gz']);
    
    %% move to permanent storage
    eval(['! mv ',proj.path.code,'tmp/sng_orient_thresh_zstatd20_',num2str(ica),'_3x3x3.nii.gz ', ... 
          proj.path.ctrl.in_ica]);

end

%% move mask and template to permanent storage
eval(['! mv ',proj.path.code,'tmp/sng_orient_mfc_mask_3x3x3.nii.gz ',proj.path.ctrl.in_ica]);
eval(['! mv ',proj.path.code,'tmp/erode2_sng_orient_mfc_mask_3x3x3.nii.gz ',proj.path.ctrl.in_ica]);
eval(['! mv ',proj.path.code,'tmp/TT_icbm452_orig.nii ',proj.path.ctrl.in_ica]);

%% clean-up
eval(['! rm ',proj.path.code,'tmp/*']);
