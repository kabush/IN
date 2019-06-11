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
logger(['Computing IN Cog Control PEL Models                  '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.ctrl.in_pel_mdl]);
    eval(['! rm -rf ',proj.path.ctrl.in_pel_mdl]);
    disp(['Creating ',proj.path.ctrl.in_pel_mdl]);
    eval(['! mkdir ',proj.path.ctrl.in_pel_mdl]);
end

%% ----------------------------------------
%% Load labels;
label_id = load([proj.path.trg.ex,'stim_ids.txt']);
v_score = load([proj.path.trg.ex,'stim_v_scores.txt']);
a_score = load([proj.path.trg.ex,'stim_a_scores.txt']);

%% Subselect intrinsic data
in_id = find(label_id==proj.param.trg.in_id);
in_v_score = v_score(in_id,1);
in_a_score = a_score(in_id,1);

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
        
        % Data is present
        data_exist = 1;
        
    catch
        logger(['   -predictions do not exist'],proj.path.logfile);
    end

    if(data_exist)
        
        % Initialize the prediction structure of this subject
        mdls = struct();

        %% ----------------------------------------
        %% For these two models, need to convert the 
        %% Predictions into probabilities first using
        p_c_v = 1./(1+exp(-zscore(in_v_score)));            %% prob of true cls
        p_c_v_all = repmat(p_c_v,1,size(prds.v_dcmp.h,2));
        p_beta_c_v = 1./(1+exp(-prds.v_dcmp.h)); %% prob of predicted cls

        %% ----------------------------------------
        %% Model Prediction Error-Likelihood
        mdls.v_dcmp.pel = abs(p_c_v_all(:,3:(end-1))-p_beta_c_v(:,3:(end-1)));
        mdls.v_indx.pel = prds.v_indx.h(:,3:(end-1));

        % save out model structure
        save([proj.path.ctrl.in_pel_mdl,subj_study,'_',name,'_mdls.mat'],'mdls');

    end
    
end
