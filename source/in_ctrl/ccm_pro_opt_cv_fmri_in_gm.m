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
logger(['Compute PRO (optional) Trajectories  '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.ctrl.in_pro_opt_mdl]);
    eval(['! rm -rf ',proj.path.ctrl.in_pro_opt_mdl]);
    disp(['Creating ',proj.path.ctrl.in_pro_opt_mdl]);
    eval(['! mkdir ',proj.path.ctrl.in_pro_opt_mdl]);
end

%% ----------------------------------------
%% ----------------------------------------
%% VALENCE inter-subject PRO predictions
affect_name = 'v';

% Predict the PRO
[pro_all_out,pro_all,mdl] = eval_pro_cv(proj,affect_name);

% Analyze the PRO
[Ntrials,Nsbjs] = size(pro_all_out);
measures = [];
predictors = [];
subjects = [];
for i=1:numel(Nsbjs)
   measures = [measures;zscore(pro_all(:,i))];
   predictors = [predictors;zscore(pro_all_out(:,i))];
   subjects = [subjects;repmat(i,Ntrials,1)];
end
tbl = table(measures,predictors,subjects,'VariableNames',{'measures','predictors','subjects'});
mdl_fe = fitlme(tbl,['measures ~ 1 + predictors']);
mdl_re = fitlme(tbl,['measures ~ 1 + predictors + (predictors|subjects)']);
fe_vs_re = compare(mdl_fe,mdl_re);
mdl = mdl_fe;
if(fe_vs_re.pValue<0.05)
    mdl=mdl_re;
    logger('  ---Random effects matters',proj.path.logfile');
end
save([proj.path.ctrl.in_pro_opt_mdl,'evaluate_pro_',affect_name,'.mat'],'mdl');



%% ----------------------------------------
%% ----------------------------------------
%% AROUSAL inter-subject PRO predictions
affect_name = 'a';

% Predict the PRO
[pro_all_out,pro_all,mdl] = eval_pro_cv(proj,affect_name);

% Analyze the PRO
[Ntrials,Nsbjs] = size(pro_all_out);
measures = [];
predictors = [];
subjects = [];
for i=1:numel(Nsbjs)
   measures = [measures;zscore(pro_all(:,i))];
   predictors = [predictors;zscore(pro_all_out(:,i))];
   subjects = [subjects;repmat(i,Ntrials,1)];
end
tbl = table(measures,predictors,subjects,'VariableNames',{'measures','predictors','subjects'});
mdl_fe = fitlme(tbl,['measures ~ 1 + predictors']);
mdl_re = fitlme(tbl,['measures ~ 1 + predictors + (predictors|subjects)']);
fe_vs_re = compare(mdl_fe,mdl_re);
mdl = mdl_fe;
if(fe_vs_re.pValue<0.05)
    mdl=mdl_re;
    logger('  ---Random effects matters',proj.path.logfile');
end
save([proj.path.ctrl.in_pro_opt_mdl,'evaluate_pro_',affect_name,'.mat'],'mdl');
