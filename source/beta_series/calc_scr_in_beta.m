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

%% Set-up Directory Structure for SCR
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.betas.scr_in_beta]);
    eval(['! rm -rf ',proj.path.betas.scr_in_beta]);
    disp(['Creating ',proj.path.betas.scr_in_beta]);
    eval(['! mkdir ',proj.path.betas.scr_in_beta]);
end

%% Load designs
load([proj.path.design,'run1_design.mat']);
load([proj.path.design,'run2_design.mat']);

%% Extract Intrinsic Stimulation Times (shifted for R5 upgrade)
run1_in_stim_times = run1_design.in_time_seq'+proj.param.trg.r5_shift;
run2_in_stim_times = run2_design.in_time_seq'+proj.param.trg.r5_shift;

%% Extract Feel Stimuluation Times
run1_feel_stim_times = [];
for i=1:numel(run1_in_stim_times)
    run1_feel_stim_times = [run1_feel_stim_times,proj.param.trg.feel_times+run1_in_stim_times(i)];
end

run2_feel_stim_times = [];
for i=1:numel(run2_in_stim_times)
    run2_feel_stim_times = [run2_feel_stim_times,proj.param.trg.feel_times+run2_in_stim_times(i)];
end

%% build design(s)
disp('bulding design files');
[prime_in_1 other_in_1] = scr_dsgn_preproc(proj,proj.param.mri.n_trs_id1,run1_in_stim_times);
[prime_in_2 other_in_2] = scr_dsgn_preproc(proj,proj.param.mri.n_trs_id2,run2_in_stim_times);
[prime_feel_1 other_feel_1] = scr_dsgn_preproc(proj,proj.param.mri.n_trs_id1,run1_feel_stim_times);
[prime_feel_2 other_feel_2] = scr_dsgn_preproc(proj,proj.param.mri.n_trs_id2,run2_feel_stim_times);

%% load subjs
subjs = load_subjs(proj);

logger(['************************************************'],proj.path.logfile);
logger(['Calculating SCR beta-series of ',num2str(numel(subjs)),' subjects'],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

grp_betas = [];

for i=1:numel(subjs)
    
    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    %% debug
    logger([subj_study,':',name],proj.path.logfile);

    %% Initialize scr beta structure
    in_betas = struct();
    feel_betas = struct();

    %% ----------------------------------------
    %% LSS of scr signal (Mumford, 2012) - Identify 1
    in_betas.id1 = [];
    feel_betas.id1 = [];
    try

        path = [proj.path.physio.scr_clean,subj_study,'_',name,'_Identify_run_1.mat'];
        load(path);

        for j=1:size(prime_in_1)
            prime = prime_in_1(j,:);
            other = other_in_1(j,:);
            mdl_in_1 = regstats(scr,[prime_in_1(j,:)',other_in_1(j,:)']);
            in_betas.id1 = [in_betas.id1,mdl_in_1.beta(2)'];
        end

        % %%Normalize
        % in_betas.id1 = zscore(in_betas.id1);

        feel_betas_tmp = [];
        for j=1:size(prime_feel_1)
            prime = prime_feel_1(j,:);
            other = other_feel_1(j,:);
            mdl_feel_1 = regstats(scr,[prime_feel_1(j,:)',other_feel_1(j,:)']);
            feel_betas_tmp = [feel_betas_tmp,mdl_feel_1.beta(2)'];
            % % old way
            % feel_betas.id1 = [feel_betas.id1,mdl_feel_1.beta(2)'];
        end
        feel_betas.id1 = reshape(feel_betas_tmp,4,15)';

    catch
        logger(['  -LSS Error: SCR of Identify run 1: ',path],proj.path.logfile);
    end

    %% ----------------------------------------
    %% LSS of scr signal (Mumford, 2012) - Identify 1
    in_betas.id2 = [];
    feel_betas.id2 = [];
    try
        
        path = [proj.path.physio.scr_clean,subj_study,'_',name,'_Identify_run_2.mat'];
        load(path);
        
        for j=1:size(prime_in_2)
            prime = prime_in_2(j,:);
            other = other_in_2(j,:);
            mdl_in_2 = regstats(scr,[prime_in_2(j,:)',other_in_2(j,:)']);
            in_betas.id2 = [in_betas.id2,mdl_in_2.beta(2)'];
        end
        
        % %%Normalize
        % in_betas.id2 = zscore(in_betas.id2);
        
        feel_betas_tmp = [];
        for j=1:size(prime_feel_2)
            prime = prime_feel_2(j,:);
            other = other_feel_2(j,:);
            mdl_feel_2 = regstats(scr,[prime_feel_2(j,:)',other_feel_2(j,:)']);
            feel_betas_tmp = [feel_betas_tmp,mdl_feel_2.beta(2)'];
            % % old way
            % feel_betas.id2 = [feel_betas.id2,mdl_feel_2.beta(2)'];
        end
        feel_betas.id2 = reshape(feel_betas_tmp,4,15)';

        % %%Normalize
        % feel_betas.id2 = zscore(feel_betas.id2);
        
    catch
        logger(['  -LSS Error: SCR of Identify run 2: ',path],proj.path.logfile);
    end

    
    %% ----------------------------------------
    %% SAVE Individual Betas
    save([proj.path.betas.scr_in_beta,subj_study,'_',name,'_in_betas.mat'],'in_betas');
    save([proj.path.betas.scr_in_beta,subj_study,'_',name,'_feel_betas.mat'],'feel_betas');
    
end
