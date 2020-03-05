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
    disp(['Removing ',proj.path.betas.emg_in_beta]);
    eval(['! rm -rf ',proj.path.betas.emg_in_beta]);
    disp(['Creating ',proj.path.betas.emg_in_beta]);
    eval(['! mkdir ',proj.path.betas.emg_in_beta]);
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

%% load subjs
subjs = load_subjs(proj);

logger(['************************************************'],proj.path.logfile);
logger(['Calculating EMG beta-series of ',num2str(numel(subjs)),' subjects'],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

for i=1:numel(subjs) % only CTM has EMG recordings that are valid
    
    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    if(strcmp(subj_study,'CTM') ~= 0 | strcmp(subj_study,'CTER') ~= 0)
        
        %% debug
        logger([subj_study,':',name],proj.path.logfile);
        
        %% Initialize emg beta structure
        in_betas = struct();
        feel_betas = struct();
        in_betas.zygo.id1 = [];
        feel_betas.zygo.id1 = [];
        in_betas.corr.id1 = [];
        feel_betas.corr.id1 = [];
        
        try
            
            path = [proj.path.physio.emg_clean,subj_study,'_',name,'_Identify_run_1_zygo.mat'];
            load(path);
            
            path = [proj.path.physio.emg_clean,subj_study,'_',name,'_Identify_run_1_corr.mat'];
            load(path);
            
            for j = 1:numel(run1_in_stim_times)
                start_time = run1_in_stim_times(j);
                start_id = start_time*proj.param.physio.hz_emg;
                stim_samples = proj.param.mri.TR*proj.param.physio.hz_emg;
                end_id = start_id+stim_samples-1;
                in_beta_zygo = sum(rect_zygo(round(start_id):round(end_id)));
                in_beta_corr = sum(rect_corr(round(start_id):round(end_id)));
                in_betas.zygo.id1 = [in_betas.zygo.id1;in_beta_zygo];
                in_betas.corr.id1 = [in_betas.corr.id1;in_beta_corr];
            end

            
            feel_beta_zygo_tmp = [];
            feel_beta_corr_tmp = [];
            for j = 1:numel(run1_feel_stim_times)
                start_time = run1_feel_stim_times(j);
                start_id = start_time*proj.param.physio.hz_emg;
                stim_samples = proj.param.mri.TR*proj.param.physio.hz_emg;
                end_id = start_id+stim_samples-1;
                feel_beta_zygo = sum(rect_zygo(round(start_id):round(end_id)));
                feel_beta_corr = sum(rect_corr(round(start_id):round(end_id)));
                feel_beta_zygo_tmp = [feel_beta_zygo_tmp;feel_beta_zygo];
                feel_beta_corr_tmp = [feel_beta_corr_tmp;feel_beta_corr];
            end
            feel_betas.zygo.id1 = reshape(feel_beta_zygo_tmp,4,15)';
            feel_betas.corr.id1 = reshape(feel_beta_corr_tmp,4,15)';
            
        catch
            logger(['  -Error: EMG of Identify run 1: ',path],proj.path.logfile);
        end
        
        %% additional emg beta structure
        in_betas.zygo.id2 = [];
        feel_betas.zygo.id2 = [];
        in_betas.corr.id2 = [];
        feel_betas.corr.id2 = [];
        
        try
            
            path = [proj.path.physio.emg_clean,subj_study,'_',name,'_Identify_run_2_zygo.mat'];
            load(path);
            
            path = [proj.path.physio.emg_clean,subj_study,'_',name,'_Identify_run_2_corr.mat'];
            load(path);
            
            for j = 1:numel(run2_in_stim_times)
                start_time = run2_in_stim_times(j);
                start_id = start_time*proj.param.physio.hz_emg;
                stim_samples = proj.param.mri.TR*proj.param.physio.hz_emg;
                end_id = start_id+stim_samples-1;
                in_beta_zygo = sum(rect_zygo(round(start_id):round(end_id)));
                in_beta_corr = sum(rect_corr(round(start_id):round(end_id)));
                in_betas.zygo.id2 = [in_betas.zygo.id2;in_beta_zygo];
                in_betas.corr.id2 = [in_betas.corr.id2;in_beta_corr];
            end
            
            feel_beta_zygo_tmp = [];
            feel_beta_corr_tmp = [];
            for j = 1:numel(run2_feel_stim_times)
                start_time = run2_feel_stim_times(j);
                start_id = start_time*proj.param.physio.hz_emg;
                stim_samples = proj.param.mri.TR*proj.param.physio.hz_emg;
                end_id = start_id+stim_samples-1;
                feel_beta_zygo = sum(rect_zygo(round(start_id):round(end_id)));
                feel_beta_corr = sum(rect_corr(round(start_id):round(end_id)));
                feel_beta_zygo_tmp = [feel_beta_zygo_tmp;feel_beta_zygo];
                feel_beta_corr_tmp = [feel_beta_corr_tmp;feel_beta_corr];
            end
            feel_betas.zygo.id2 = reshape(feel_beta_zygo_tmp,4,15)';
            feel_betas.corr.id2 = reshape(feel_beta_corr_tmp,4,15)';
            
        catch
            logger(['  -Error: EMG of Identify run 2: ',path],proj.path.logfile);
        end
        
        %% ----------------------------------------
        %% SAVE Individual Betas
        save([proj.path.betas.emg_in_beta,subj_study,'_',name,'_in_betas.mat'],'in_betas');
        save([proj.path.betas.emg_in_beta,subj_study,'_',name,'_feel_betas.mat'],'feel_betas');
        
    end

end
