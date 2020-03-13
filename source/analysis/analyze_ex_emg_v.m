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
    disp(['Removing ',proj.path.analysis.ex_emg_v]);
    eval(['! rm -rf ',proj.path.analysis.ex_emg_v]);
    disp(['Creating ',proj.path.analysis.ex_emg_v]);
    eval(['! mkdir ',proj.path.analysis.ex_emg_v]);
end

logger(['*************************************************'],proj.path.logfile);
logger(['Analyzing IN (Zygomaticus)          '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
calc_ex_emg_v(proj,'zygo');

logger(['*************************************************'],proj.path.logfile);
logger(['Analyzing IN (Corrugator)          '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
calc_ex_emg_v(proj,'corr');
