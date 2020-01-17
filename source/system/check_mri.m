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
logger(' MRI Quality Check ');
logger('----------------------------------------');

%% Create the subjects to be analyzed (possible multiple studies)
subjs = proj.process.subjs; 

%% Preprocess fMRI of each subject in subjects list 
for i=1:numel(subjs)
    
    subj = subjs{i};
    
    %% extract subject info
    subj_study = subj.study;
    name = subj.name;
    
    % files are usable
    subj.mri.ok = 0;
    subj.mri.anat.ok = 0;
    subj.mri.id1.ok = 0;
    subj.mri.id2.ok = 0;
    subj.mri.gm_mask.ok = 0;
    
    % exist
    subj.mri.anat.exist = 0;
    subj.mri.id1.exist = 0;
    subj.mri.id2.exist = 0;
    
    % motion
    subj.mri.id1.censor = [];
    subj.mri.id2.censor = [];
    subj.mri.id1.fd = [];
    subj.mri.id2.fd = [];
    subj.mri.id1.censor_ok = 0;
    subj.mri.id2.censor_ok = 0;
    subj.mri.id1.motion_ok = 0;
    subj.mri.id2.motion_ok = 0;
    
    % check for existence of anat file
    if exist([proj.path.mri.mri_clean,subj_study,'_',name,'/anat/' ...
              ,subj_study,'.',name,'.anat.seg.fsl.MNI.GM+tlrc.BRIK'],'file')==2
        subj.mri.anat.exist = 1;
    end
    
    % check for existence of mri file (identify 1)
    if exist([proj.path.mri.mri_clean,subj_study,'_',name,'/identify/run1/' ...
              ,subj_study,'.',name,'.identify.run1.scaled.resid.nii'],'file')==2
        subj.mri.id1.exist = 1;
        
        % get censor and fd
        if(subj.mri.id1.exist)
            subj.mri.id1.censor = load([proj.path.mri.mri_clean,subj_study,'_',name,'/identify/run1/' ...
                                ,subj_study,'.',name,'.identify.run1.censor.1D']);
            subj.mri.id1.fd = load([proj.path.mri.mri_clean,subj_study,'_',name,'/identify/run1/' ...
                                ,subj_study,'.',name,'.identify.run1.FD.1D']);

            % check censor file length
            if(numel(subj.mri.id1.censor)==proj.param.mri.n_trs_id1)
                subj.mri.id1.censor_ok = 1;
            end
            
            % calculate motion
            N_all = numel(subj.mri.id1.fd);
            N_bad = numel(find(subj.mri.id1.fd>proj.param.mri.FD_thresh));
            if((N_bad/N_all)<proj.param.mri.FD_bad_frac)
                subj.mri.id1.motion_ok=1;
            end
            
            %debug
            logger(['Identify task motion (fd>0.5): ',num2str(round(100*N_bad/N_all)),'%'],proj.path.logfile);
            
        end
        
    end
    
    % check for existence of mri file (identify 2)
    if exist([proj.path.mri.mri_clean,subj_study,'_',name,'/identify/run2/' ...
              ,subj_study,'.',name,'.identify.run2.scaled.resid.nii'],'file')==2
        subj.mri.id2.exist = 1;
        
        % get censor and fd
        if(subj.mri.id2.exist)
            subj.mri.id2.censor = load([proj.path.mri.mri_clean,subj_study,'_',name,'/identify/run2/' ...
                                ,subj_study,'.',name,'.identify.run2.censor.1D']);
            subj.mri.id2.fd = load([proj.path.mri.mri_clean,subj_study,'_',name,'/identify/run2/' ...
                                ,subj_study,'.',name,'.identify.run2.FD.1D']);
            

            % check censor file length
            if(numel(subj.mri.id2.censor)==proj.param.mri.n_trs_id2)
                subj.mri.id2.censor_ok = 1;
            end
            
            % calculate motion
            N_all = numel(subj.mri.id2.fd);
            N_bad = numel(find(subj.mri.id2.fd>proj.param.mri.FD_thresh));
            if(N_bad/N_all<proj.param.mri.FD_bad_frac)
                subj.mri.id2.motion_ok=1;
            end
            
        end
        
    end
    
    % check for existence of mri mask
    if exist([proj.path.mri.gm_mask,subj_study,'.',name,'.gm.nii'],'file')==2
        subj.mri.gm_mask.ok = 1;
    end
    
    %% ----------------------------------------
    %% Quality control logic
    
    %% Anat images
    if(subj.mri.anat.exist)
        subj.mri.anat.ok=1;
    end
    
    %% EPI images
    if(subj.mri.id1.exist & subj.mri.id1.censor_ok & subj.mri.id1.motion_ok)
        subj.mri.id1.ok = 1;
    end
    
    if(subj.mri.id2.exist & subj.mri.id2.censor_ok & subj.mri.id2.motion_ok)
        subj.mri.id2.ok = 1;
    end
    
    %% Combined anat and EPI check for subject
    % 1) check that all fmri components are available (4 pieces)
    if(subj.mri.anat.ok+subj.mri.id1.ok+subj.mri.id2.ok+subj.mri.gm_mask.ok==4)
        subj.mri.ok = 1;
    end
    
    % assign master flag to false
    if(subj.mri.ok==0)
        subj.ok = 0;
        logger([' -Excluding from study: ',subj_study,'_',name],proj.path.logfile);
    end
    
    % output
    proj.process.subjs{i} = subj;
    
end

%% Indicate quality check has been completed
proj.check.mri = 1;

%% ----------------------------------------
%% Write out amended project structure
save('proj.mat','proj');