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

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.betas.fmri_in_beta]);
    eval(['! rm -rf ',proj.path.betas.fmri_in_beta]);
    disp(['Creating ',proj.path.betas.fmri_in_beta]);
    eval(['! mkdir ',proj.path.betas.fmri_in_beta]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

logger(['************************************************'],proj.path.logfile);
logger(['Calculating fMRI beta-series (IN) of ',num2str(numel(subjs)),' subjects'],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

for i=1:numel(subjs)
    
    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    task = 'identify';  %%% HARDCODED FOR THIS STEP

    %% log subject
    logger([subj_study,':',name],proj.path.logfile);

    %% ----------------------------------------
    %% Load and concatenate censor files (save to tmp)
    censor1 = load([proj.path.mri.mri_clean,subj_study,'_',name,'/',task,'/run1/',subj_study,'.',name,'.',task,'.run1.censor.1D']);
    censor2 = load([proj.path.mri.mri_clean,subj_study,'_',name,'/',task,'/run2/',subj_study,'.',name,'.',task,'.run2.censor.1D']);
    censor = [censor1;censor2];
    save([proj.path.code,'tmp/',subj_study,'_',name,'.',task,'.cmb.censor.1D'],'censor','-ascii');

    %% ----------------------------------------
    %% Load and concatenate motion regression files (Power et al.,
    %% 20??) (save to tmp)
    motion1 = load([proj.path.mri.mri_clean,subj_study,'_',name,'/',task,'/run1/',subj_study,'.',name,'.',task,'.run1.motion.1D']);
    motion2 = load([proj.path.mri.mri_clean,subj_study,'_',name,'/',task,'/run2/',subj_study,'.',name,'.',task,'.run2.motion.1D']);
    motion = [motion1;motion2];
    save([proj.path.code,'tmp/',subj_study,'_',name,'.',task,'.cmb.motion.1D'],'motion','-ascii');

    motion_square1 = load([proj.path.mri.mri_clean,subj_study,'_',name,'/',task,'/run1/',subj_study,'.',name,'.',task,'.run1.motion.square.1D']);
    motion_square2 = load([proj.path.mri.mri_clean,subj_study,'_',name,'/',task,'/run2/',subj_study,'.',name,'.',task,'.run2.motion.square.1D']);
    motion_square = [motion_square1;motion_square2];
    save([proj.path.code,'tmp/',subj_study,'_',name,'.',task,'.cmb.motion.square.1D'],'motion_square','-ascii');
    
    motion_pre_t1 = load([proj.path.mri.mri_clean,subj_study,'_',name,'/',task,'/run1/',subj_study,'.',name,'.',task,'.run1.motion_pre_t.1D']);
    motion_pre_t2 = load([proj.path.mri.mri_clean,subj_study,'_',name,'/',task,'/run2/',subj_study,'.',name,'.',task,'.run2.motion_pre_t.1D']);
    motion_pre_t = [motion_pre_t1;motion_pre_t2];
    save([proj.path.code,'tmp/',subj_study,'_',name,'.',task,'.cmb.motion_pre_t.1D'],'motion_pre_t','-ascii');
    
    motion_pre_t_square1 = load([proj.path.mri.mri_clean,subj_study,'_',name,'/',task,'/run1/',subj_study,'.',name,'.',task,'.run1.motion_pre_t_square.1D']);
    motion_pre_t_square2 = load([proj.path.mri.mri_clean,subj_study,'_',name,'/',task,'/run2/',subj_study,'.',name,'.',task,'.run2.motion_pre_t_square.1D']);
    motion_pre_t_square = [motion_pre_t_square1;motion_pre_t_square2];
    save([proj.path.code,'tmp/',subj_study,'_',name,'.',task,'.cmb.motion_pre_t_square.1D'],'motion_pre_t_square','-ascii');
    
    %% ----------------------------------------
    %% Write important project data to tmp for use by afni code
    dlmwrite([proj.path.code,'tmp/target_path.txt'],proj.path.trg.in,'');
    dlmwrite([proj.path.code,'tmp/project_name.txt'],proj.path.name,'');
    dlmwrite([proj.path.code,'tmp/out_path.txt'],proj.path.betas.fmri_in_beta,'');
    dlmwrite([proj.path.code,'tmp/task.txt'],task,'');
    dlmwrite([proj.path.code,'tmp/TR.txt'],proj.param.mri.TR,'');
    dlmwrite([proj.path.code,'tmp/nTR.txt'],proj.param.mri.n_trs_id1,'');
    dlmwrite([proj.path.code,'tmp/stim_t.txt'],proj.param.trg.stim_t,'');
    dlmwrite([proj.path.code,'tmp/study.txt'],subj_study,'');
    dlmwrite([proj.path.code,'tmp/subject.txt'],name,'');

    %% ----------------------------------------
    %% fit the beta-series
    eval(['! ',proj.path.code,'source/beta_series/fmri_3dlss ',proj.path.home,' ',proj.path.name]);
    
    %% ----------------------------------------
    %% Clean-up
    eval(['! rm ',proj.path.code,'tmp/*']);
    
end
