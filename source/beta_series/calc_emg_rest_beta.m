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
    disp(['Removing ',proj.path.betas.emg_rest_beta]);
    eval(['! rm -rf ',proj.path.betas.emg_rest_beta]);
    disp(['Creating ',proj.path.betas.emg_rest_beta]);
    eval(['! mkdir ',proj.path.betas.emg_rest_beta]);
end

%% Locally rename project params
N_trs = proj.param.mri.n_trs_rest;
N_sample = 30; 
N_trans = proj.param.rest.n_trs_trans;
N_tail = proj.param.rest.n_trs_tail;

%% Extract Sim. IN Stimuluation Times
rest_in_stim_times = proj.param.mri.TR*randsample((N_trans+1):(N_trs-N_tail),N_sample);

%% Extract Sim. Feel Stimuluation Times
rest_feel_stim_times = [];
for i=1:numel(rest_in_stim_times)
    rest_feel_stim_times = [rest_feel_stim_times,proj.param.trg.feel_times+rest_in_stim_times(i)];
end

%% load subjs
subjs = load_subjs(proj);

logger(['************************************************'],proj.path.logfile);
logger(['Calculating EMG REST beta-series of ',num2str(numel(subjs)),' subjects'],proj.path.logfile);
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
        in_betas.zygo.rest = [];
        feel_betas.zygo.rest = [];
        in_betas.corr.rest = [];
        feel_betas.corr.rest = [];
        
        try
            
            path = [proj.path.physio.emg_clean,subj_study,'_',name,'_Rest_zygo.mat'];
            load(path);
            
            path = [proj.path.physio.emg_clean,subj_study,'_',name,'_Rest_corr.mat'];
            load(path);
            
            for j = 1:numel(rest_in_stim_times)
                start_time = rest_in_stim_times(j);
                start_id = start_time*proj.param.physio.hz_emg;
                stim_samples = proj.param.mri.TR*proj.param.physio.hz_emg;
                end_id = start_id+stim_samples-1;
                in_beta_zygo = sum(rect_zygo(round(start_id):round(end_id)));
                in_beta_corr = sum(rect_corr(round(start_id):round(end_id)));
                in_betas.zygo.rest = [in_betas.zygo.rest;in_beta_zygo];
                in_betas.corr.rest = [in_betas.corr.rest;in_beta_corr];
            end

            feel_beta_zygo_tmp = [];
            feel_beta_corr_tmp = [];
            for j = 1:numel(rest_feel_stim_times)
                start_time = rest_feel_stim_times(j);
                start_id = start_time*proj.param.physio.hz_emg;
                stim_samples = proj.param.mri.TR*proj.param.physio.hz_emg;
                end_id = start_id+stim_samples-1;
                feel_beta_zygo = sum(rect_zygo(round(start_id):round(end_id)));
                feel_beta_corr = sum(rect_corr(round(start_id):round(end_id)));
                feel_beta_zygo_tmp = [feel_beta_zygo_tmp;feel_beta_zygo];
                feel_beta_corr_tmp = [feel_beta_corr_tmp;feel_beta_corr];
            end
            feel_betas.zygo.rest = reshape(feel_beta_zygo_tmp,4,30)';
            feel_betas.corr.rest = reshape(feel_beta_corr_tmp,4,30)';
            
        catch
            logger(['  -Error: EMG of Rest: ',path],proj.path.logfile);
        end
        
        %% ----------------------------------------
        %% SAVE Individual Betas
        save([proj.path.betas.emg_rest_beta,subj_study,'_',name,'_in_betas.mat'],'in_betas');
        save([proj.path.betas.emg_rest_beta,subj_study,'_',name,'_feel_betas.mat'],'feel_betas');
        
    end

end
