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
    disp(['Removing ',proj.path.betas.scr_rest_beta]);
    eval(['! rm -rf ',proj.path.betas.scr_rest_beta]);
    disp(['Creating ',proj.path.betas.scr_rest_beta]);
    eval(['! mkdir ',proj.path.betas.scr_rest_beta]);
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

%% build design(s)
disp('bulding design files');
[prime_in other_in] = scr_dsgn_preproc(proj,proj.param.mri.n_trs_rest,rest_in_stim_times);
[prime_feel other_feel] = scr_dsgn_preproc(proj,proj.param.mri.n_trs_rest,rest_feel_stim_times);

%% load subjs
subjs = load_subjs(proj);

logger(['************************************************'],proj.path.logfile);
logger(['Calculating SCR REST beta-series of ',num2str(numel(subjs)),' subjects'],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

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
    %% LSS of scr signal (Mumford, 2012) - REST
    in_betas = [];
    feel_betas = [];
    try

        path = [proj.path.physio.scr_clean,subj_study,'_',name,'_Rest.mat'];
        load(path);

        for j=1:size(prime_in)
            prime = prime_in(j,:);
            other = other_in(j,:);
            mdl_in = regstats(scr,[prime_in(j,:)',other_in(j,:)']);
            in_betas = [in_betas,mdl_in.beta(2)'];
        end
        
        feel_betas_tmp = [];
        for j=1:size(prime_feel)
            prime = prime_feel(j,:);
            other = other_feel(j,:);
            mdl_feel = regstats(scr,[prime_feel(j,:)',other_feel(j,:)']);
            feel_betas_tmp = [feel_betas_tmp,mdl_feel.beta(2)'];
        end
        feel_betas = reshape(feel_betas_tmp,4,30)';  


    catch
        logger(['  -LSS Error: SCR of Rest: ',path],proj.path.logfile);
    end

    %% ----------------------------------------
    %% SAVE Individual Betas
    save([proj.path.betas.scr_rest_beta,subj_study,'_',name,'_in_betas.mat'],'in_betas');
    save([proj.path.betas.scr_rest_beta,subj_study,'_',name,'_feel_betas.mat'],'feel_betas');
    

end
