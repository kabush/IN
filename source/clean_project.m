%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

%% ========================================
%% Load in path data
load('proj.mat');

if(proj.flag.clean_build)

    %% Set-up base data pathway
    disp(['Removing ',proj.path.data]);
    eval(['! rm -rf ',proj.path.data]);

    %% Create project directory
    disp(['Creating ',proj.path.data,' and all sub-directories']);
    eval(['! mkdir ',proj.path.data]);

    %% Create all top-level directories
    eval(['! mkdir ',proj.path.data,proj.path.system.name]);
    eval(['! mkdir ',proj.path.data,proj.path.analysis.name]);
    eval(['! mkdir ',proj.path.data,proj.path.betas.name]);
    eval(['! mkdir ',proj.path.data,proj.path.ctrl.name]);
    eval(['! mkdir ',proj.path.data,proj.path.haufe.name]);
    eval(['! mkdir ',proj.path.data,proj.path.mri.name]);
    eval(['! mkdir ',proj.path.data,proj.path.mvpa.name]);
    eval(['! mkdir ',proj.path.data,proj.path.physio.name]);
    eval(['! mkdir ',proj.path.data,proj.path.trg.name]);

end