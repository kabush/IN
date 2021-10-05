%%========================================
%%========================================
%%
%% Keith Bush, PhD (2019)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

function [act_mu,act_std,err_mu,err_std] = compute_dsc_act_err_stats(proj,affect_name,act_5part)

%% Initialize log section
logger(['*************************************************'],proj.path.logfile);
if(strcmp(affect_name,'v'))
    logger(['Computing Discrete Actions (VALENCE)'],proj.path.logfile);
else
    logger(['Computing Discrete Actions (AROUSAL)'],proj.path.logfile);
end
logger(['*************************************************'],proj.path.logfile);


subjs = load_subjs(proj);

dsc_act_all = [];
dsc_err_all = [];

for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    % log processing of subject
    logger([subj_study,'_',name],proj.path.logfile);

    data_exist = 0;
    try

        %% load dynamics
        load([proj.path.ctrl.in_dyn,subj_study,'_',name,'_prds.mat']);
        
        %% Data is present
        data_exist = 1;
        
    catch
        logger(['   -predictions do not exist'],proj.path.logfile);
    end

    if(data_exist)

        %% ******* ACTIONS ********
        post_vals = eval(['prds.',affect_name,'_dcmp.h(:,4:end)']);
        pre_vals = eval(['prds.',affect_name,'_dcmp.h(:,3:(end-1))']);
        actions = post_vals-pre_vals;
        actions_1d = zscore(reshape(actions',1,prod(size(actions))));

        dsc_actions_1d = 0*actions_1d;
        dsc_actions_1d(find(actions_1d>act_5part(3)))=1; 
        dsc_actions_1d(find(actions_1d>act_5part(4)))=2; 
        dsc_actions_1d(find(actions_1d<act_5part(2)))=-1; 
        dsc_actions_1d(find(actions_1d<act_5part(1)))=-2; 
        dsc_act_all = [dsc_act_all,dsc_actions_1d];

        %% ******* ERRORS ********
        errors = (prds.v_dcmp.err(:,4:end));
        errors_1d = zscore(reshape(errors',1,prod(size(errors))));
        dsc_err_all = [dsc_err_all,errors_1d];
        
    end

end

% Compute group params
act_mu = mean(dsc_act_all);
act_std = std(dsc_act_all);
err_mu = mean(dsc_err_all);
err_std = std(dsc_err_all);
