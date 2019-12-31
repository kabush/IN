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

% Observe Q-function parameters 
[q_perf_v,act_err_v,sig_test_v] = calc_q_param_perf(proj,...
                                                  discount_set,...
                                                  reward_frac_set,...
                                                  Q_traj_all,...
                                                  Q_rand_all,...
                                                  act_err_all);

%% ----------------------------------------
%% AROUSAL analysis

% Load Q-function performance
load([proj.path.ctrl.in_evc_opt_mdl,'Q_traj_all_a.mat']);
load([proj.path.ctrl.in_evc_opt_mdl,'Q_rand_all_a.mat']);
load([proj.path.ctrl.in_evc_opt_mdl,'act_err_all_a.mat']);

% Observe Q-function parameters 
[q_perf_a,act_err_a,sig_test_a] = calc_q_param_perf(proj,...
                                                  discount_set,...
                                                  reward_frac_set,...
                                                  Q_traj_all,...
                                                  Q_rand_all,...
                                                  act_err_all);

