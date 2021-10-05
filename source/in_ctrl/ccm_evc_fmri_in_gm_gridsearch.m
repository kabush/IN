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
logger(['Gridsearch of IN EVC Parameters  '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.ctrl.in_evc_opt_mdl]);
    eval(['! rm -rf ',proj.path.ctrl.in_evc_opt_mdl]);
    disp(['Creating ',proj.path.ctrl.in_evc_opt_mdl]);
    eval(['! mkdir ',proj.path.ctrl.in_evc_opt_mdl]);
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

%%%% RUN ALL SUBJECTS 
Nrand = numel(subjs);
rand_subj_ids = 1:numel(subjs);  

% Meta RL Parameter search
discount_set = proj.param.ctrl.discount_set; % raw discount value used in Q-function update
reward_frac_set = proj.param.ctrl.reward_frac_set; % balance between reward/action

% Calculate set-sizes
Ndsct = numel(discount_set);
Nfrac = numel(reward_frac_set);

%% ----------------------------------------
%% ----------------------------------------
%% VALENCE Computations for a random subject set

logger(['*************************************************'],proj.path.logfile);
logger(['Gridsearch EVC Parameters (VALENCE) '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

Q_traj_all = zeros(Ndsct,Nfrac,Nrand,30,4);
Q_rand_all = zeros(Ndsct,Nfrac,Nrand,30,4);
act_err_all = zeros(Ndsct,Nfrac,Nrand,30,4);

grp_act_mu = act_mu_v;
grp_act_std = act_std_v;
grp_err_mu = err_mu_v;
grp_err_std = err_std_v;
act_5part = act_5part_v;
affect_name = 'v';

for b=1:Nfrac

    for a=1:Ndsct


        gamma = discount_set(a);
        rwrd_act_f = reward_frac_set(b);

        logger(['gamma: ',num2str(gamma),', frac: ',num2str(rwrd_act_f)],proj.path.logfile);

        [Q_traj_cv,Q_rand_cv,act_err_cv,eval_subj_ids] = ...
            eval_qfunc_param_cv(proj, ...
                                act_5part, ...
                                grp_act_mu, ...
                                grp_act_std, ...
                                grp_err_mu, ...
                                grp_err_std, ...
                                rand_subj_ids,gamma, ...
                                rwrd_act_f,affect_name);

        Q_traj_all(a,b,:,:,:) = Q_traj_cv;
        Q_rand_all(a,b,:,:,:) = Q_rand_cv;
        act_err_all(a,b,:,:,:) = act_err_cv;

        % Save intermediate results
        save([proj.path.ctrl.in_evc_opt_mdl,'Q_traj_all_',affect_name,'.mat'],'Q_traj_all');
        save([proj.path.ctrl.in_evc_opt_mdl,'Q_rand_all_',affect_name,'.mat'],'Q_rand_all');
        save([proj.path.ctrl.in_evc_opt_mdl,'act_err_all_',affect_name,'.mat'],'act_err_all');
        
        % Save subjects involved in calculations
        save([proj.path.ctrl.in_evc_opt_mdl,'rand_subj_ids_',affect_name,'.mat'],'rand_subj_ids');
        save([proj.path.ctrl.in_evc_opt_mdl,'eval_subj_ids_',affect_name,'.mat'],'eval_subj_ids');

        % Remove the output of the fitted Q-iteration
        eval(['! rm ',proj.path.ctrl.in_evc_opt_mdl,'*_result_',affect_name,'.mat']);
        
    end

end

%% ----------------------------------------
%% ----------------------------------------
%% AROUSAL Computations for a random subject set

logger(['*************************************************'],proj.path.logfile);
logger(['Gridsearch EVC Parameters (AROUSAL) '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

Q_traj_all = zeros(Ndsct,Nfrac,Nrand,30,4);
Q_rand_all = zeros(Ndsct,Nfrac,Nrand,30,4);
act_err_all = zeros(Ndsct,Nfrac,Nrand,30,4);

grp_act_mu = act_mu_a;
grp_act_std = act_std_a;
grp_err_mu = err_mu_a;
grp_err_std = err_std_a;
act_5part = act_5part_a;
affect_name = 'a';

for b=1:Nfrac

    for a=1:Ndsct

        gamma = discount_set(a);
        rwrd_act_f = reward_frac_set(b);

        logger(['gamma: ',num2str(gamma),', frac: ',num2str(rwrd_act_f)],proj.path.logfile);

        [Q_traj_cv,Q_rand_cv,act_err_cv,eval_subj_ids] = ...
            eval_qfunc_param_cv(proj, ...
                                act_5part, ...
                                grp_act_mu, ...
                                grp_act_std, ...
                                grp_err_mu, ...
                                grp_err_std, ...
                                rand_subj_ids,gamma, ...
                                rwrd_act_f,affect_name);

        Q_traj_all(a,b,:,:,:) = Q_traj_cv;
        Q_rand_all(a,b,:,:,:) = Q_rand_cv;
        act_err_all(a,b,:,:,:) = act_err_cv;

        % Save intermediate results
        save([proj.path.ctrl.in_evc_opt_mdl,'Q_traj_all_',affect_name,'.mat'],'Q_traj_all');
        save([proj.path.ctrl.in_evc_opt_mdl,'Q_rand_all_',affect_name,'.mat'],'Q_rand_all');
        save([proj.path.ctrl.in_evc_opt_mdl,'act_err_all_',affect_name,'.mat'],'act_err_all');
        
        % Save subjects involved in calculations
        save([proj.path.ctrl.in_evc_opt_mdl,'rand_subj_ids_',affect_name,'.mat'],'rand_subj_ids');
        save([proj.path.ctrl.in_evc_opt_mdl,'eval_subj_ids_',affect_name,'.mat'],'eval_subj_ids');

        % Remove the output of the fitted Q-iteration
        eval(['! rm ',proj.path.ctrl.in_evc_opt_mdl,'*_result_',affect_name,'.mat']);
        
    end

end
