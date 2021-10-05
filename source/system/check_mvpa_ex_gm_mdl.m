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
logger(' MVPA EX GM Model Quality Check ');
logger('----------------------------------------');

%% Create the subjects to be analyzed (possible multiple studies)
subjs = proj.process.subjs;

if(proj.process.mvpa_ex_gm_mdl)

    %% Preprocess fMRI of each subject in subjects list 
    for i=1:numel(subjs)
       
        subj = subjs{i};
        
        %% extract subject info
        subj_study = subj.study;
        name = subj.name;
        
        % files are usable quality
        subj.mvpa.ex_gm_mdl.ok = 0;

        subj.mvpa.ex_gm_mdl_v.ok = 0;
        subj.mvpa.ex_gm_mdl_a.ok = 0;

        subj.mvpa.ex_gm_mdl_v.exist = 0;
        subj.mvpa.ex_gm_mdl_a.exist = 0;

        subj.mvpa.ex_gm_mdl_v.nan_ok = 0;
        subj.mvpa.ex_gm_mdl_a.nan_ok = 0;        

        % ----------------------------------------
        % check for existence of model (Valence)
        if exist([proj.path.mvpa.fmri_ex_gm_mdl,subj_study,'_',name,'_v_model.mat']);
            subj.mvpa.ex_gm_mdl_v.exist = 1;
        end

        % check for nans
        if(subj.mvpa.ex_gm_mdl_v.exist);

            % load model
            load([proj.path.mvpa.fmri_ex_gm_mdl,subj_study,'_',name,'_v_model.mat']);

            % check valence predictions for NAN
             if(numel(find(isnan(v_model.Beta))))
                 disp([subj_study,'_',name,': ', ...
                       num2str(numel(find(isnan(v_model.Beta)))),' VAL ' ...
                       'model betas) NAN']);
            else
                subj.mvpa.ex_gm_mdl_v.nan_ok = 1;
            end

        end


        % ----------------------------------------
        % check for existence of model (Arousal)
        if exist([proj.path.mvpa.fmri_ex_gm_mdl,subj_study,'_',name,'_a_model.mat']);
            subj.mvpa.ex_gm_mdl_a.exist = 1;
        end

        % check for nans
        if(subj.mvpa.ex_gm_mdl_a.exist);

            % load model
            load([proj.path.mvpa.fmri_ex_gm_mdl,subj_study,'_',name,'_a_model.mat']);

            % check valence predictions for NAN
             if(numel(find(isnan(a_model.Beta))))
                 disp([subj_study,'_',name,': ', ...
                       num2str(numel(find(isnan(a_model.Beta)))),' ARO ' ...
                       'model betas) NAN']);
            else
                subj.mvpa.ex_gm_mdl_a.nan_ok = 1;
            end

        end
       
        %% ----------------------------------------
        %% Quality control logic

        %% valence model
        if(subj.mvpa.ex_gm_mdl_v.exist & subj.mvpa.ex_gm_mdl_v.nan_ok)
            subj.mvpa.ex_gm_mdl_v.ok = 1;
        end

        %% arousal model
        if(subj.mvpa.ex_gm_mdl_a.exist & subj.mvpa.ex_gm_mdl_a.nan_ok)
            subj.mvpa.ex_gm_mdl_a.ok = 1;
        end

        %% total prediction
        if(subj.mvpa.ex_gm_mdl_v.ok & subj.mvpa.ex_gm_mdl_a.ok)
            subj.mvpa.ex_gm_mdl.ok = 1;
        end

        % assign master flag to false
        if(subj.mvpa.ex_gm_mdl.ok==0)
            subj.ok = 0;
        end
        
        % output
        proj.process.subjs{i} = subj;
        
    end

else
    disp('No MVPA EX GS CLS predictions to quality check.  Processing not complete.');
end

%% Indicate quality check has been completed
proj.check.mvpa_ex_gm_mdl = 1;

%% ----------------------------------------
%% Write out amended project structure
save('proj.mat','proj');