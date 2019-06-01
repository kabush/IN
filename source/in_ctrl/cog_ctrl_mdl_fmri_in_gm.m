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

%% Initialize log section
logger(['*************************************************'],proj.path.logfile);
logger(['Computing IN Cog Control Models                  '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.ctrl.in_ccm]);
    eval(['! rm -rf ',proj.path.ctrl.in_ccm]);
    disp(['Creating ',proj.path.ctrl.in_ccm]);
    eval(['! mkdir ',proj.path.ctrl.in_ccm]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% ----------------------------------------
%% Transform beta-series into affect series {v,a}
for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    % log processing of subject
    logger([subj_study,'_',name],proj.path.logfile);

    data_exist = 0;
    try

        % load dynamics
        load([proj.path.ctrl.in_dyn,subj_study,'_',name,'_prds.mat']);
        
        %% Data is present
        data_exist = 1;
        
    catch
        logger(['   -predictions do not exist'],proj.path.logfile);
    end

    if(data_exist)

        %% Initialize the prediction structure of this subject
        mdls = struct();

        %% ********************
        %% *** DONE TO HERE ***
        %% ********************

        %% ----------------------------------------
        %%  1) Error
        mdls.v_dcmp.err = prds.v_dcmp.err;
        mdls.v_indx.err = prds.v_indx.err;  %% error here, not indices

        %% ----------------------------------------
        %%  2) Conflict


        %% ----------------------------------------
        %% For these two models, need to convert the 
        %% Predictions into probabilities first using
        %% Platt scaling
        %%  3) Prediction Error-Likelihood
        %%  4) PRO


        %% ----------------------------------------
        %% For this model to work we need to get the random
        %% forests model of RL predictions
        %%  5) EVC

    end
    
end
