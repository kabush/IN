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
logger([' Analyzing VR Q values (GRID SEARCH)             '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% Meta RL Parameter search
discount_set = [0:.1:1];
reward_frac_set = [0:.2:1]; % balance between reward/action

%% ----------------------------------------
%% VALENCE analysis

% Load Q-function performance
load([proj.path.ctrl.in_evc_opt_mdl,'Q_traj_all_v.mat']);
load([proj.path.ctrl.in_evc_opt_mdl,'Q_rand_all_v.mat']);
load([proj.path.ctrl.in_evc_opt_mdl,'act_err_all_v.mat']);

% Observe Q-function parameters impacts on control
[q_perf,act_err,sig_test] = calc_q_param_perf(proj,...
                                              discount_set,...
                                              reward_frac_set,...
                                              Q_traj_all,...
                                              Q_rand_all,...
                                              act_err_all);

% Solve for optimal paramter
[gamma,frac] = calc_q_param_opt(discount_set,reward_frac_set,act_err);
gamma
frac

% Save out findings
save([proj.path.ctrl.in_evc_opt_mdl,'q_perf_v.mat'],'q_perf');
save([proj.path.ctrl.in_evc_opt_mdl,'act_err_v.mat'],'act_err');
save([proj.path.ctrl.in_evc_opt_mdl,'sig_test_v.mat'],'sig_test');
save([proj.path.ctrl.in_evc_opt_mdl,'gamma_v.mat'],'gamma');
save([proj.path.ctrl.in_evc_opt_mdl,'frac_v.mat'],'frac');

%% ----------------------------------------
%% AROUSAL analysis

% Load Q-function performance
load([proj.path.ctrl.in_evc_opt_mdl,'Q_traj_all_a.mat']);
load([proj.path.ctrl.in_evc_opt_mdl,'Q_rand_all_a.mat']);
load([proj.path.ctrl.in_evc_opt_mdl,'act_err_all_a.mat']);

% Observe Q-function parameters impacts on control
[q_perf,act_err,sig_test] = calc_q_param_perf(proj,...
                                              discount_set,...
                                              reward_frac_set,...
                                              Q_traj_all,...
                                              Q_rand_all,...
                                              act_err_all);

[gamma,frac] = calc_q_param_opt(discount_set,reward_frac_set,act_err);
gamma
frac

% Save out findings
save([proj.path.ctrl.in_evc_opt_mdl,'q_perf_a.mat'],'q_perf');
save([proj.path.ctrl.in_evc_opt_mdl,'act_err_a.mat'],'act_err');
save([proj.path.ctrl.in_evc_opt_mdl,'sig_test_a.mat'],'sig_test');
save([proj.path.ctrl.in_evc_opt_mdl,'gamma_a.mat'],'gamma');
save([proj.path.ctrl.in_evc_opt_mdl,'frac_a.mat'],'frac');