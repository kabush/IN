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
    disp(['Removing ',proj.path.analysis.emg_rest_entrain]);
    eval(['! rm -rf ',proj.path.analysis.emg_rest_entrain]);
    disp(['Creating ',proj.path.analysis.emg_rest_entrain]);
    eval(['! mkdir ',proj.path.analysis.emg_rest_entrain]);
end

%% Initialize log section
logger(['*************************************************'],proj.path.logfile);
logger(['Analyzing REST Entrain (Zygomaticus)          '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
calc_rest_emg_entrain(proj,'zygo');

logger(['*************************************************'],proj.path.logfile);
logger(['Analyzing REST Entrain (Corrugator)          '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
calc_rest_emg_entrain(proj,'corr');
