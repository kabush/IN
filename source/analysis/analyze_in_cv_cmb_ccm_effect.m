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
    disp(['Removing ',proj.path.analysis.in_cv_cmb_ccm_effect]);
    eval(['! rm -rf ',proj.path.analysis.in_cv_cmb_ccm_effect]);
    disp(['Creating ',proj.path.analysis.in_cv_cmb_ccm_effect]);
    eval(['! mkdir ',proj.path.analysis.in_cv_cmb_ccm_effect]);
end

logger(['*************************************************'],proj.path.logfile);
logger(['Estimate Relative CCM Effects (VALENCE)     '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
calc_in_cv_cmb_ccm_effect(proj,'v');

logger(['*************************************************'],proj.path.logfile);
logger(['Estimate Relative CCM Effects (AROUSAL)     '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
calc_in_cv_cmb_ccm_effect(proj,'a');

logger(['*************************************************'],proj.path.logfile);
logger(['Compare CCM Effects across affect     '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
cmp_in_cv_cmb_ccm_effect(proj);

