%%========================================
%%========================================
%%
%% Keith Bush, PhD (2019)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

%% Load in path data
load('proj.mat');

%% Initialize log section
logger(['*************************************************'],proj.path.logfile);
logger(['Compute EVC Q-values (via Cross-validation)  '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.ctrl.in_evc_cv_mdl]);
    eval(['! rm -rf ',proj.path.ctrl.in_evc_cv_mdl]);
    disp(['Creating ',proj.path.ctrl.in_evc_cv_mdl]);
    eval(['! mkdir ',proj.path.ctrl.in_evc_cv_mdl]);
end

%% ----------------------------------------
%% Compute group partitions of 5-bin discrete action space
[act_5part_v] = compute_action_partition(proj,'v');
[act_5part_a] = compute_action_partition(proj,'a');

%% ----------------------------------------
%% Compute group discrete action and error stats
[act_mu_v,act_std_v,err_mu_v,err_std_v] = compute_dsc_act_err_stats(proj,'v',act_5part_v);
[act_mu_a,act_std_a,err_mu_a,err_std_a] = compute_dsc_act_err_stats(proj,'a',act_5part_a);

subjs = load_subjs(proj);

for i=1:numel(subjs)

    % extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;
    
    % log processing of subject
    logger([subj_study,'_',name],proj.path.logfile);
    
    mdls = struct();

    %% ----------------------------------------
    %% VALENCE Computations for a random subject set
    
    logger(['Compute EVC (VALENCE) '],proj.path.logfile);
    grp_act_mu = act_mu_v;
    grp_act_std = act_std_v;
    grp_err_mu = err_mu_v;
    grp_err_std = err_std_v;
    act_5part = act_5part_v;
    affect_name = 'v';
    
    % Load in optimal parameters
    load([proj.path.ctrl.in_evc_opt_mdl,'gamma_v.mat']); 
    load([proj.path.ctrl.in_evc_opt_mdl,'frac_v.mat']); 
    
    % Fit subjects Q-functions (& trajs) with optimal params
    [dcmp,indx] = eval_qfunc_cv(proj, ...
                                subjs{i}, ...
                                act_5part, ...
                                grp_act_mu, ...
                                grp_act_std, ...
                                grp_err_mu, ...
                                grp_err_std, ...
                                gamma, ...
                                frac,...
                                affect_name);
    mdls.v_dcmp = dcmp;
    mdls.v_indx.evc = indx;
    
    %% ----------------------------------------
    %% AROUSAL Computations for a random subject set
    
    logger(['Compute EVC (AROUSAL) '],proj.path.logfile);
    grp_act_mu = act_mu_a;
    grp_act_std = act_std_a;
    grp_err_mu = err_mu_a;
    grp_err_std = err_std_a;
    act_5part = act_5part_a;
    affect_name = 'a';
    
    % Load in optimal parameters
    load([proj.path.ctrl.in_evc_opt_mdl,'gamma_a.mat']); 
    load([proj.path.ctrl.in_evc_opt_mdl,'frac_a.mat']);
    
    % Fit subjects Q-functions (& trajs) with optimal params
    [dcmp,indx] = eval_qfunc_cv(proj, ...
                                subjs{i}, ...
                                act_5part, ...
                                grp_act_mu, ...
                                grp_act_std, ...
                                grp_err_mu, ...
                                grp_err_std, ...
                                gamma, ...
                                frac, ...
                                affect_name);
    mdls.a_dcmp = dcmp;
    mdls.a_indx = indx;
    
    %% ----------------------------------------
    %% save models
    save([proj.path.ctrl.in_evc_cv_mdl,subj_study,'_',name,'_mdls.mat'],'mdls');
    
end
