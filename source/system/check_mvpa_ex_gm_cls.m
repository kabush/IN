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

logger('----------------------------------------');
logger(' MVPA EX GM CLS Quality Check ');
logger('----------------------------------------');

%% Create the subjects to be analyzed (possible multiple studies)
subjs = proj.process.subjs;

%% Preprocess fMRI of each subject in subjects list 
for i=1:numel(subjs)
   
    subj = subjs{i};
    
    %% extract subject info
    subj_study = subj.study;
    name = subj.name;
    
    % files are usable quality
    subj.mvpa.ex_gm_cls.ok = 0;

    subj.mvpa.ex_gm_cls_v.ok = 0;
    subj.mvpa.ex_gm_cls_a.ok = 0;

    subj.mvpa.ex_gm_cls.exist = 0;
    subj.mvpa.ex_gm_cls_v.nan_ok = 0;
    subj.mvpa.ex_gm_cls_a.nan_ok = 0;        

    % check for existence of EX GM CLS predictions
    if exist([proj.path.mvpa.fmri_ex_gm_cls,subj_study,'_',name,'_prds.mat']);
        subj.mvpa.ex_gm_cls.exist = 1;
    end

    % check for nans
    if(subj.mvpa.ex_gm_cls.exist);

        % load predictions
        load([proj.path.mvpa.fmri_ex_gm_cls,subj_study,'_',name,'_prds.mat']);

        % check valence predictions for NAN
        prds_mu = mean(prds.v_cls_acc,1);
        if(numel(find(isnan(prds_mu))))
            disp([subj_study,'_',name,': ', ...
                  num2str(numel(find(isnan(prds_mu)))),' VAL ' ...
                                'predictions(s) NAN']);
        else
            subj.mvpa.ex_gm_cls_v.nan_ok = 1;
        end

        % check arousal predictions for NAN
        prds_mu = mean(prds.a_cls_acc,1);
        if(numel(find(isnan(prds_mu))))
            disp([subj_study,'_',name,': ', ...
                  num2str(numel(find(isnan(prds_mu)))),' ARO ' ...
                                'predictions(s) NAN']);
        else
            subj.mvpa.ex_gm_cls_a.nan_ok = 1;
        end

    end
   
    %% ----------------------------------------
    %% Quality control logic

    %% predictions valence
    if(subj.mvpa.ex_gm_cls.exist & subj.mvpa.ex_gm_cls_v.nan_ok)
        subj.mvpa.ex_gm_cls_v.ok = 1;
    end

    %% predictions arousal
    if(subj.mvpa.ex_gm_cls.exist & subj.mvpa.ex_gm_cls_a.nan_ok)
        subj.mvpa.ex_gm_cls_a.ok = 1;
    end

    %% total prediction
    if(subj.mvpa.ex_gm_cls_v.ok & subj.mvpa.ex_gm_cls_a.ok)
        subj.mvpa.ex_gm_cls.ok = 1;
    end

    % assign master flag to false
    if(subj.mvpa.ex_gm_cls.ok==0)
        subj.ok = 0;
    end
    
    % output
    proj.process.subjs{i} = subj;
    
end

%% Indicate quality check has been completed
proj.check.mvpa_ex_gm_cls = 1;

%% ----------------------------------------
%% Write out amended project structure
save('proj.mat','proj');