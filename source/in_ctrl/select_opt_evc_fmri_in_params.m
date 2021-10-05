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

discount_set = proj.param.ctrl.discount_set; 
reward_frac_set = proj.param.ctrl.reward_frac_set; 

%% ----------------------------------------
%% VALENCE analysis

% Load Q-function performance
load([proj.path.ctrl.in_evc_opt_mdl,'Q_traj_all_v.mat']);
load([proj.path.ctrl.in_evc_opt_mdl,'Q_rand_all_v.mat']);
load([proj.path.ctrl.in_evc_opt_mdl,'act_err_all_v.mat']);

% Load subjects
subjs = load_subjs(proj);

% Find usable entries in data structures (b/c of poorly
% written (and slow) output from the gridsearch
usable_cnt = 0;
for i=1:numel(subjs)
    mu = mean(mean(squeeze(Q_traj_all(1,1,i,:,:))));
    if(abs(mu)>0)
        usable_cnt = i;
    end
end
Nsbj = usable_cnt;

%% conduct out-of-sample parameter estimation
all_gamma = [];
all_frac = [];
for i=1:Nsbj

    subj_ids = setdiff(1:Nsbj,i);
    
    % Observe Q-function parameters impacts on control
    [q_perf,act_err,sig_test] = calc_q_param_perf(proj,...
                                                  discount_set,...
                                                  reward_frac_set,...
                                                  Q_traj_all(:,:,subj_ids,:,:),...
                                                  Q_rand_all(:,:,subj_ids,:,:),...
                                                  act_err_all(:,:,subj_ids,:,:));
    
    % Solve for optimal paramter
    [gamma,frac] = calc_q_param_opt(discount_set,reward_frac_set,act_err);
    disp(['gamma=',num2str(gamma),', frac=',num2str(frac)]);            
    
    all_gamma = [all_gamma,gamma];
    all_frac = [all_frac,frac];

end

%% find closest parameter to cross-validated group mean
gamma_diff = sqrt((discount_set - mean(all_gamma)).^2);
gamma = discount_set(find(gamma_diff==min(gamma_diff)));
frac_diff = sqrt((reward_frac_set - mean(all_frac)).^2);
frac = reward_frac_set(find(frac_diff==min(frac_diff)));
disp(['GRP CV MU: gamma=',num2str(gamma),', frac=',num2str(frac)]);            

% Observe Q-function parameters impacts on control
[q_perf,act_err,sig_test] = calc_q_param_perf(proj,...
                                              discount_set,...
                                              reward_frac_set,...
                                              Q_traj_all(:,:,1:Nsbj,:,:),...
                                              Q_rand_all(:,:,1:Nsbj,:,:),...
                                              act_err_all(:,:,1:Nsbj,:,:));

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

% Load subjects
subjs = load_subjs(proj);

% Find usable entries in data structures (b/c of poorly
% written (and slow) output from the gridsearch
usable_cnt = 0;
for i=1:numel(subjs)
    mu = mean(mean(squeeze(Q_traj_all(1,1,i,:,:))));
    if(abs(mu)>0)
        usable_cnt = i;
    end
end
Nsbj = usable_cnt;

%% conduct out-of-sample parameter estimation
all_gamma = [];
all_frac = [];
for i=1:Nsbj

    subj_ids = setdiff(1:Nsbj,i);
    
    % Observe Q-function parameters impacts on control
    [q_perf,act_err,sig_test] = calc_q_param_perf(proj,...
                                                  discount_set,...
                                                  reward_frac_set,...
                                                  Q_traj_all(:,:,subj_ids,:,:),...
                                                  Q_rand_all(:,:,subj_ids,:,:),...
                                                  act_err_all(:,:,subj_ids,:,:));
    
    % Solve for optimal paramter
    [gamma,frac] = calc_q_param_opt(discount_set,reward_frac_set,act_err);
    disp(['gamma=',num2str(gamma),', frac=',num2str(frac)]);            
    
    all_gamma = [all_gamma,gamma];
    all_frac = [all_frac,frac];

end

%% find closest parameter to cross-validated group mean
gamma_diff = sqrt((discount_set - mean(all_gamma)).^2);
gamma = discount_set(find(gamma_diff==min(gamma_diff)));
frac_diff = sqrt((reward_frac_set - mean(all_frac)).^2);
frac = reward_frac_set(find(frac_diff==min(frac_diff)));
disp(['GRP CV MU: gamma=',num2str(gamma),', frac=',num2str(frac)]);            

% Observe Q-function parameters impacts on control
[q_perf,act_err,sig_test] = calc_q_param_perf(proj,...
                                              discount_set,...
                                              reward_frac_set,...
                                              Q_traj_all(:,:,1:Nsbj,:,:),...
                                              Q_rand_all(:,:,1:Nsbj,:,:),...
                                              act_err_all(:,:,1:Nsbj,:,:));

% Save out findings
save([proj.path.ctrl.in_evc_opt_mdl,'q_perf_a.mat'],'q_perf');
save([proj.path.ctrl.in_evc_opt_mdl,'act_err_a.mat'],'act_err');
save([proj.path.ctrl.in_evc_opt_mdl,'sig_test_a.mat'],'sig_test');
save([proj.path.ctrl.in_evc_opt_mdl,'gamma_a.mat'],'gamma');
save([proj.path.ctrl.in_evc_opt_mdl,'frac_a.mat'],'frac');
