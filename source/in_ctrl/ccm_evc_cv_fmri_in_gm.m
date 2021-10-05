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
logger(['Compute EVC Trajectories  '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.ctrl.in_evc_icv_mdl]);
    eval(['! rm -rf ',proj.path.ctrl.in_evc_icv_mdl]);
    disp(['Creating ',proj.path.ctrl.in_evc_icv_mdl]);
    eval(['! mkdir ',proj.path.ctrl.in_evc_icv_mdl]);
end

%% ----------------------------------------
%% Compute group partitions of 5-bin discrete action space
[act_5part_v] = compute_action_partition(proj,'v');
[act_5part_a] = compute_action_partition(proj,'a');

%% ----------------------------------------
%% Compute group discrete action and error stats
[act_mu_v,act_std_v,err_mu_v,err_std_v] = compute_dsc_act_err_stats(proj,'v',act_5part_v);
[act_mu_a,act_std_a,err_mu_a,err_std_a] = compute_dsc_act_err_stats(proj,'a',act_5part_a);

%% ----------------------------------------
%% Compute Q-functions for all random subject set
%% to characterize EVC params.

% load subjs
subjs = load_subjs(proj);
Nsubj = numel(subjs);
subj_ids = 1:Nsubj;

%% ----------------------------------------
%% ----------------------------------------
%% VALENCE intersubject predictions

logger(['*************************************************'],proj.path.logfile);
logger(['Compute inter-subject EVC Parameters (VALENCE) '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

grp_act_mu = act_mu_v;
grp_act_std = act_std_v;
grp_err_mu = err_mu_v;
grp_err_std = err_std_v;
act_5part = act_5part_v;
affect_name = 'v';

% Load in optimal parameters
load([proj.path.ctrl.in_evc_opt_mdl,'gamma_v.mat']); 
load([proj.path.ctrl.in_evc_opt_mdl,'frac_v.mat']); 

logger(['gamma: ',num2str(gamma),', frac: ',num2str(frac)],proj.path.logfile);

% Run inter-subject cross-validation
[Q_cv,Q_traj_cv,Q_rand_cv,Q_bst_cv,Q_cnf_cv,act_err_cv,] = ...
    eval_qfunc_cv(proj, ...
                  act_5part, ...
                  grp_act_mu, ...
                  grp_act_std, ...
                  grp_err_mu, ...
                  grp_err_std, ...
                  subj_ids,gamma, ...
                  frac,affect_name);

% Save intermediate results
save([proj.path.ctrl.in_evc_icv_mdl,'Q_cv_',affect_name,'.mat'],'Q_cv');
save([proj.path.ctrl.in_evc_icv_mdl,'Q_traj_',affect_name,'.mat'],'Q_traj_cv');
save([proj.path.ctrl.in_evc_icv_mdl,'Q_rand_',affect_name,'.mat'],'Q_rand_cv');
save([proj.path.ctrl.in_evc_icv_mdl,'Q_bst_',affect_name,'.mat'],'Q_bst_cv');
save([proj.path.ctrl.in_evc_icv_mdl,'Q_cnf_',affect_name,'.mat'],'Q_cnf_cv');
save([proj.path.ctrl.in_evc_icv_mdl,'act_err_',affect_name,'.mat'],'act_err_cv');

%% ----------------------------------------
%% ----------------------------------------
%% AROUSAL intersubject predictions

logger(['*************************************************'],proj.path.logfile);
logger(['Compute inter-subject EVC Parameters (AROUSAL) '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

grp_act_mu = act_mu_a;
grp_act_std = act_std_a;
grp_err_mu = err_mu_a;
grp_err_std = err_std_a;
act_5part = act_5part_a;
affect_name = 'a';

% Load in optimal parameters
load([proj.path.ctrl.in_evc_opt_mdl,'gamma_a.mat']); 
load([proj.path.ctrl.in_evc_opt_mdl,'frac_a.mat']); 

logger(['gamma: ',num2str(gamma),', frac: ',num2str(frac)],proj.path.logfile);

% Run inter-subject cross-validation
[Q_cv,Q_traj_cv,Q_rand_cv,Q_bst_cv,Q_cnf_cv,act_err_cv,] = ...
    eval_qfunc_cv(proj, ...
                  act_5part, ...
                  grp_act_mu, ...
                  grp_act_std, ...
                  grp_err_mu, ...
                  grp_err_std, ...
                  subj_ids,gamma, ...
                  frac,affect_name);

% Save intermediate results
save([proj.path.ctrl.in_evc_icv_mdl,'Q_cv_',affect_name,'.mat'],'Q_cv');
save([proj.path.ctrl.in_evc_icv_mdl,'Q_traj_',affect_name,'.mat'],'Q_traj_cv');
save([proj.path.ctrl.in_evc_icv_mdl,'Q_rand_',affect_name,'.mat'],'Q_rand_cv');
save([proj.path.ctrl.in_evc_icv_mdl,'Q_bst_',affect_name,'.mat'],'Q_bst_cv');
save([proj.path.ctrl.in_evc_icv_mdl,'Q_cnf_',affect_name,'.mat'],'Q_cnf_cv');
save([proj.path.ctrl.in_evc_icv_mdl,'act_err_',affect_name,'.mat'],'act_err_cv');
