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
    disp(['Removing ',proj.path.analysis.in_cv_cmb_3dlme]);
    eval(['! rm -rf ',proj.path.analysis.in_cv_cmb_3dlme]);
    disp(['Creating ',proj.path.analysis.in_cv_cmb_3dlme]);
    eval(['! mkdir ',proj.path.analysis.in_cv_cmb_3dlme]);
end

%% Initialize log section
logger(['*************************************************'],proj.path.logfile);
logger(['3dLME detection of Cog Ctrl Voxels (VALENCE)     '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
calc_in_cv_cmb_fmri_3dlme(proj,'v');

logger(['*************************************************'],proj.path.logfile);
logger(['3dLME detection of Cog Ctrl Voxels (AROUSAL)     '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
calc_in_cv_cmb_fmri_3dlme(proj,'a');
