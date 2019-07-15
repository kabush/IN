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
logger(' MRI Beta-series Quality Check ');
logger('----------------------------------------');

%% Create the subjects to be analyzed (possible multiple studies)
subjs = proj.process.subjs; %load_subjs(proj);

%% Preprocess fMRI of each subject in subjects list 
for i=1:numel(subjs)  %nans @ 2, 18, 32
   
    subj = subjs{i};
    
    %% extract subject info
    subj_study = subj.study;
    name = subj.name;
    
    % files are usable quality
    subj.beta.mri_ex_id.ok = 0; %usable
    subj.beta.mri_ex_id.exist = 0; %exist
    subj.beta.mri_ex_id.nan_ok = 0; %nan corruption
    subj.beta.mri_ex_id.nan_ids = [];
    
    % check for existence of Identify beta-series
    if exist([proj.path.betas.fmri_ex_beta,subj_study,'_',name,'_lss.nii']);
        subj.beta.mri_ex_id.exist = 1;
    end

    % check for nans
    if(subj.beta.mri_ex_id.exist)

        base_nii = load_nii([proj.path.betas.fmri_ex_beta,subj_study,'_',name,'_lss.nii']);
        base_img = vec_img_2d_nii(base_nii);
       
        base_mu = mean(base_img);
        if(numel(find(isnan(base_mu))))

            disp([subj_study,'_',name,': ', ...
                  num2str(numel(find(isnan(base_mu)))),' volume(s) NAN']);
            subj.beta.mri_ex_id.nan_ids = find(isnan(base_mu));

        else
            subj.beta.mri_ex_id.nan_ok = 1;
        end
    end
   
    %% ----------------------------------------
    %% Quality control logic

    %% beta-series
    if(subj.beta.mri_ex_id.exist & subj.beta.mri_ex_id.nan_ok)
        subj.beta.mri_ex_id.ok = 1;
    end

    % assign master flag to false
    if(subj.beta.mri_ex_id.ok==0)
        subj.ok = 0;
    end
    
    % output
    proj.process.subjs{i} = subj;
    
end

%% Indicate quality check has been completed
proj.check.beta_mri_ex_id = 1;

%% ----------------------------------------
%% Write out amended project structure
save('proj.mat','proj');