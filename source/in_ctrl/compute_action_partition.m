%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

%% ----------------------------------------
%% discretize action set

function [act_5part] = compute_action_partition(proj,affect_name)

%% Initialize log section
logger(['*************************************************'],proj.path.logfile);
if(strcmp(affect_name,'v'))
    logger(['Computing Action Partition (VALENCE)'],proj.path.logfile);
else
    logger(['Computing Action Partition (AROUSAL)'],proj.path.logfile);
end
logger(['*************************************************'],proj.path.logfile);

%% load subjs
subjs = load_subjs(proj);

Nsbj = 0;
act_all = [];
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
        
        %% Data is present
        data_exist = 1;
        
    catch
        logger(['   -predictions do not exist'],proj.path.logfile);
    end

    if(data_exist)

        post_vals = eval(['prds.',affect_name,'_dcmp.h(:,4:end)']);
        pre_vals = eval(['prds.',affect_name,'_dcmp.h(:,3:(end-1))']);
        actions = post_vals-pre_vals;
        actions_1d = zscore(reshape(actions',1,prod(size(actions))));       
        act_all = [act_all,actions_1d];

        %% count subjects
        Nsbj = Nsbj + 1;

    end

end

%% compute standard deviation
act_std = std(act_all);

%% compute partition based on standard deviations
act_5part = [-2*act_std,-act_std,act_std,2*act_std];