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

%% ----------------------------------------
%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.analysis.in_cv_cmb_clust_thresh]);
    eval(['! rm -rf ',proj.path.analysis.in_cv_cmb_clust_thresh]);
    disp(['Creating ',proj.path.analysis.in_cv_cmb_clust_thresh]);
    eval(['! mkdir ',proj.path.analysis.in_cv_cmb_clust_thresh]);
end

%% ----------------------------------------
%% VALENCE ANATOMY
%% ----------------------------------------

% estimate cluster sizes
logger(['*************************************************'],proj.path.logfile);
logger(['3dLME cluster threshold estimation (VALENCE)     '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
calc_in_cv_cmb_clust_thresh(proj,'v');

% z-stat
logger(['*************************************************'],proj.path.logfile);
logger(['Apply 3dLME cluster threshold (VALENCE)   '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
apply_in_cv_cmb_clust_thresh(proj,'v');

logger(['*************************************************'],proj.path.logfile);
logger(['Apply mFC Mask (VALENCE)   '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
apply_in_cv_cmb_mfc_mask(proj,'v');

logger(['*************************************************'],proj.path.logfile);
logger(['Apply dACC Mask (VALENCE)   '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
apply_in_cv_cmb_dacc_mask(proj,'v');

% f-stat
logger(['*************************************************'],proj.path.logfile);
logger(['Apply 3dLME cluster threshold (VALENCE)   '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
apply_in_cv_cmb_clust_thresh_2fstat(proj,'v');

logger(['*************************************************'],proj.path.logfile);
logger(['Apply mFC Mask to FSTAT (VALENCE)'],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
apply_in_cv_cmb_mfc_mask_2fstat(proj,'v');

logger(['*************************************************'],proj.path.logfile);
logger(['Apply dACC Mask to FSTAT (VALENCE)'],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
apply_in_cv_cmb_dacc_mask_2fstat(proj,'v');

%% ----------------------------------------
%% AROUSAL ANATOMY
%% ----------------------------------------

% estimate cluster sizes
logger(['*************************************************'],proj.path.logfile);
logger(['3dLME cluster threshold estimation (AROUSAL)     '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
calc_in_cv_cmb_clust_thresh(proj,'a');

% z-stat
logger(['*************************************************'],proj.path.logfile);
logger(['Apply 3dLME cluster threshold (AROUSAL)   '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
apply_in_cv_cmb_clust_thresh(proj,'a');

logger(['*************************************************'],proj.path.logfile);
logger(['Apply mFC Mask (AROUSAL)   '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
apply_in_cv_cmb_mfc_mask(proj,'a');

logger(['*************************************************'],proj.path.logfile);
logger(['Apply dACC Mask (AROUSAL)   '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
apply_in_cv_cmb_dacc_mask(proj,'a');

% f-stat
logger(['*************************************************'],proj.path.logfile);
logger(['Apply 3dLME cluster threshold (AROUSAL)   '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
apply_in_cv_cmb_clust_thresh_2fstat(proj,'a');

logger(['*************************************************'],proj.path.logfile);
logger(['Apply mFC Mask to FSTAT (AROUSAL)'],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
apply_in_cv_cmb_mfc_mask_2fstat(proj,'a');

logger(['*************************************************'],proj.path.logfile);
logger(['Apply dACC Mask to FSTAT (AROUSAL)'],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
apply_in_cv_cmb_dacc_mask_2fstat(proj,'a');
