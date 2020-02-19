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
    disp(['Removing ',proj.path.analysis.in_base_clust_thresh]);
    eval(['! rm -rf ',proj.path.analysis.in_base_clust_thresh]);
    disp(['Creating ',proj.path.analysis.in_base_clust_thresh]);
    eval(['! mkdir ',proj.path.analysis.in_base_clust_thresh]);
end

logger(['*************************************************'],proj.path.logfile);
logger(['3dLME cluster threshold estimation (VALENCE)     '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
calc_in_base_clust_thresh(proj,'v');

logger(['*************************************************'],proj.path.logfile);
logger(['Apply 3dLME cluster threshold z-tests (VALENCE)   '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
apply_in_base_clust_thresh(proj,'v');

logger(['*************************************************'],proj.path.logfile);
logger(['Apply mFC Mask z-tests (VALENCE)   '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
apply_in_base_mfc_mask(proj,'v');

%% ****************************************

logger(['*************************************************'],proj.path.logfile);
logger(['3dLME cluster threshold estimation (AROUSAL)     '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
calc_in_base_clust_thresh(proj,'a');

logger(['*************************************************'],proj.path.logfile);
logger(['Apply 3dLME cluster threshold (AROUSAL)  '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
apply_in_base_clust_thresh(proj,'a');
 
logger(['*************************************************'],proj.path.logfile);
logger(['Apply mFC Mask (AROUSAL)   '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
apply_in_base_mfc_mask(proj,'a');

