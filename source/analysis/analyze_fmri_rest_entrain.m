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

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.analysis.fmri_rest_entrain]);
    eval(['! rm -rf ',proj.path.analysis.fmri_rest_entrain]);
    disp(['Creating ',proj.path.analysis.fmri_rest_entrain]);
    eval(['! mkdir ',proj.path.analysis.fmri_rest_entrain]);
end

%% Initialize log section
logger(['*************************************************'],proj.path.logfile);
logger(['Analyzing REST Entrain (VALENCE)          '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
calc_fmri_rest_entrain(proj,'v');

logger(['*************************************************'],proj.path.logfile);
logger(['Analyzing REST Entrain (AROUSAL)          '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
calc_fmri_rest_entrain(proj,'a');
