%%========================================
%%========================================
%%
%% Keith Bush, PhD (2017)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

%% Load in path data
load('proj.mat');

%% Initialize log section
logger(['********************************************'],proj.path.logfile);
logger([' Analyzing GS Classification (Refit)        '],proj.path.logfile);
logger(['********************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.analysis.gs_cls_refit]);
    eval(['! rm -rf ',proj.path.analysis.gs_cls_refit]);
    disp(['Creating ',proj.path.analysis.gs_cls_refit]);
    eval(['! mkdir ',proj.path.analysis.gs_cls_refit]);
end

% %% ----------------------------------------
% %% Load labels;
% v_label = load([proj.path.trg.ex,'stim_v_labs.txt']);
% a_label = load([proj.path.trg.ex,'stim_a_labs.txt']);
% label_id = load([proj.path.trg.ex,'stim_ids.txt']);
% v_score = load([proj.path.trg.ex,'stim_v_scores.txt']);
% a_score = load([proj.path.trg.ex,'stim_a_scores.txt']);

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);



all_raw_acc_v = [];
all_raw_acc_a = [];

for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;
    
    disp([subj_study,'_',name]);
    
    % Load classification performance
    load([proj.path.mvpa.fmri_ex_gs_cls,subj_study,'_',name,'_prds.mat']);
    
    all_raw_acc_v = [all_raw_acc_v;mean(prds.v_cls_acc)];
    all_raw_acc_a = [all_raw_acc_a;mean(prds.a_cls_acc)];

end

%% ----------------------------------------
%% "Knows-what-it-knows" performance

%% ----------------------------------------
%% VALENCE Analysis
grp_v_corr_perf = [];
grp_v_incorr_perf = [];

all_corr_ids_v = {};
all_incorr_ids_v = {};
all_refit_ids_v = {};

for i=1:numel(subjs)

    %% Test sub-selection based on group performance
    val_ids = setdiff(1:numel(subjs),i);
    
    %% Test group performance against prior
    corr_ids_v = [];
    incorr_ids_v = [];
    refit_ids_v = [];

    for j=1:size(all_raw_acc_v,2)
        
        %% Compute fraction correct
        img_frac = mean(all_raw_acc_v(val_ids,j));
        
        %% ----------------------------------------
        %% Binomial test (corr/incorrect)
        [phat ci] = binofit((numel(subjs)-1)/2,numel(subjs)-1,0.05);
        
        %% Store correct/incorrect ids
        if(img_frac>ci(2))
            corr_ids_v = [corr_ids_v,j];
        end
        if(img_frac<ci(1))
            incorr_ids_v = [incorr_ids_v,j];
        end
        
        %% ----------------------------------------
        %% Binomial test refit
        [phat ci] = binofit(50,100,0.05);
        if(img_frac>ci(2))
            refit_ids_v = [refit_ids_v,j];
        end
        
        
    end
    
    %% Store all ids
    all_corr_ids_v{i} = corr_ids_v;
    all_incorr_ids_v{i} = incorr_ids_v;
    all_refit_ids_v{i} = refit_ids_v;
    
    %% Compute adjusted subject classification performance
    grp_v_corr_perf = [grp_v_corr_perf;mean(all_raw_acc_v(i,corr_ids_v))];
    grp_v_incorr_perf = [grp_v_incorr_perf;mean(all_raw_acc_v(i,incorr_ids_v))];
    
end    

%% report findings for this thresh level
[h p ci stats] = ttest(grp_v_corr_perf,0.5);
[numel(corr_ids_v),mean(mean(all_raw_acc_v,2)),mean(grp_v_corr_perf),p]


%% ----------------------------------------
%% AROUSAL Analysis
grp_a_corr_perf = [];
grp_a_incorr_perf = [];

all_corr_ids_a = {};
all_incorr_ids_a = {};
all_refit_ids_a = {};

for i=1:numel(subjs)

    %% Test sub-selection based on group performance
    val_ids = setdiff(1:numel(subjs),i);
    
    %% Test group performance against prior
    corr_ids_a = [];
    incorr_ids_a = [];
    refit_ids_a = [];

    for j=1:size(all_raw_acc_a,2)
        
        %% Compute fraction correct
        img_frac = mean(all_raw_acc_a(val_ids,j));
        
        %% ----------------------------------------
        %% Binomial test (corr/incorrect)
        [phat ci] = binofit((numel(subjs)-1)/2,numel(subjs)-1,0.05);
        
        %% Store correct/incorrect ids
        if(img_frac>ci(2))
            corr_ids_a = [corr_ids_a,j];
        end
        if(img_frac<ci(1))
            incorr_ids_a = [incorr_ids_a,j];
        end
        
        %% ----------------------------------------
        %% Binomial test refit
        [phat ci] = binofit(50,100,0.05);
        if(img_frac>ci(2))
            refit_ids_a = [refit_ids_a,j];
        end
        
        
    end
    
    %% Store all ids
    all_corr_ids_a{i} = corr_ids_a;
    all_incorr_ids_a{i} = incorr_ids_a;
    all_refit_ids_a{i} = refit_ids_a;
    
    %% Compute adjusted subject classification performance
    grp_a_corr_perf = [grp_a_corr_perf;mean(all_raw_acc_a(i,corr_ids_a))];
    grp_a_incorr_perf = [grp_a_incorr_perf;mean(all_raw_acc_a(i,incorr_ids_a))];
    
end    


%% ----------------------------------------
%% Valence output

disp(['Grp VAL accuracy, pre-refit=',num2str(mean(mean(all_raw_acc_v,2))),...
      ', post-refit=',num2str(mean(grp_v_corr_perf))]);
      
ncorr = numel(corr_ids_v);
[a b] = binofit(ncorr/2,ncorr,0.05);
sscnt_v = numel(find(grp_v_corr_perf>b(2)));

disp(['Sing. Subj VAL significant post-refit=',num2str(sscnt_v),'/',num2str(numel(subjs))]);

%% ----------------------------------------
%% Valence output
disp(['Grp ARO accuracy, pre-refit=',num2str(mean(mean(all_raw_acc_a,2))),...
      ', post-refit=',num2str(mean(grp_a_corr_perf))]);
      
ncorr = numel(corr_ids_a);
[a b] = binofit(ncorr/2,ncorr,0.05);
sscnt_a = numel(find(grp_a_corr_perf>b(2)));

disp(['Sing. Subj ARO significant post-refit=',num2str(sscnt_a),'/',num2str(numel(subjs))]);
