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
    disp(['Removing ',proj.path.analysis.in_base_3dlme]);
    eval(['! rm -rf ',proj.path.analysis.in_base_3dlme]);
    disp(['Creating ',proj.path.analysis.in_base_3dlme]);
    eval(['! mkdir ',proj.path.analysis.in_base_3dlme]);
end

%% Initialize log section
logger(['*************************************************'],proj.path.logfile);
logger(['3dLME detection of base activation (VALENCE)     '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
calc_in_base_fmri_3dlme(proj,'v');

logger(['*************************************************'],proj.path.logfile);
logger(['3dLME detection of base activation (AROUSAL)      '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
calc_in_base_fmri_3dlme(proj,'a');
