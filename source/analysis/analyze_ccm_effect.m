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
    disp(['Removing ',proj.path.analysis.in_ccm_effect]);
    eval(['! rm -rf ',proj.path.analysis.in_ccm_effect]);
    disp(['Creating ',proj.path.analysis.in_ccm_effect]);
    eval(['! mkdir ',proj.path.analysis.in_ccm_effect]);
end

%% Initialize log section
logger(['*************************************************'],proj.path.logfile);
logger(['Compute CCM Prediction Effects (VALENCE)       '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
calc_ccm_effect(proj,'v');

% %% Initialize log section
% logger(['*************************************************'],proj.path.logfile);
% logger(['Compute CCM Prediction Effects (AROUSAL)       '],proj.path.logfile);
% logger(['*************************************************'],proj.path.logfile);
% calc_ccm_effect(proj,'a');
% 
% %% Group comparison here
% % (TBD)


