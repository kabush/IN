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
    disp(['Removing ',proj.path.analysis.vr_skill]);
    eval(['! rm -rf ',proj.path.analysis.vr_skill]);
    disp(['Creating ',proj.path.analysis.vr_skill]);
    eval(['! mkdir ',proj.path.analysis.vr_skill]);
end

%% Initialize log section
logger(['*************************************************'],proj.path.logfile);
logger(['Analyzing ER Skill (VALENCE)          '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
calc_in_skill(proj,'v');

logger(['*************************************************'],proj.path.logfile);
logger(['Analyzing ER Skill (AROUSAL)          '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
calc_in_skill(proj,'a');
