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
logger(['Computing IN Cog Control Error Models            '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.ctrl.in_err_mdl]);
    eval(['! rm -rf ',proj.path.ctrl.in_err_mdl]);
    disp(['Creating ',proj.path.ctrl.in_err_mdl]);
    eval(['! mkdir ',proj.path.ctrl.in_err_mdl]);
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
        
        % data is present
        data_exist = 1;
        
    catch
        logger(['   -predictions do not exist'],proj.path.logfile);
    end

    if(data_exist)

        % Initialize the prediction structure of this subject
        mdls = struct();

        %% ----------------------------------------
        %% Model Error
        
        %valence
        mdls.v_dcmp = prds.v_dcmp.err(:,3:(end-1));
        mdls.v_indx = prds.v_indx.err(:,3:(end-1)); 

        %arousal
        mdls.a_dcmp = prds.a_dcmp.err(:,3:(end-1));
        mdls.a_indx = prds.a_indx.err(:,3:(end-1)); 

        % save out model structure
        save([proj.path.ctrl.in_err_mdl,subj_study,'_',name,'_mdls.mat'],'mdls');

    end
    
end
