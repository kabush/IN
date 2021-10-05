%%========================================
%%========================================
%%
%% Keith Bush, PhD (2020)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

%% Load in path data
load('proj.mat');

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.analysis.in_emg]);
    eval(['! rm -rf ',proj.path.analysis.in_emg]);
    disp(['Creating ',proj.path.analysis.in_emg]);
    eval(['! mkdir ',proj.path.analysis.in_emg]);
end

%% Initialize log section
logger(['*************************************************'],proj.path.logfile);
logger(['Analyzing IN (Zygomaticus)          '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
calc_in_emg(proj,'zygo');

logger(['*************************************************'],proj.path.logfile);
logger(['Analyzing IN (Corrugator)          '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
calc_in_emg(proj,'corr');
