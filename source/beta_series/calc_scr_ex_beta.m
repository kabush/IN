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
    disp(['Removing ',proj.path.betas.scr_ex_beta]);
    eval(['! rm -rf ',proj.path.betas.scr_ex_beta]);
    disp(['Creating ',proj.path.betas.scr_ex_beta]);
    eval(['! mkdir ',proj.path.betas.scr_ex_beta]);
end

%% Load designs
load([proj.path.design,'run1_design.mat']);
load([proj.path.design,'run2_design.mat']);

%% Extract Extrinsic Stimulation Times (shifted for R5 upgrade)
run1_ex_stim_times = run1_design.ex_time_seq'+proj.param.trg.r5_shift; 
run2_ex_stim_times = run2_design.ex_time_seq'+proj.param.trg.r5_shift;

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
[prime_ex_1 other_ex_1] = scr_dsgn_preproc(proj,proj.param.mri.n_trs_id1,run1_ex_stim_times);
[prime_ex_2 other_ex_2] = scr_dsgn_preproc(proj,proj.param.mri.n_trs_id2,run2_ex_stim_times);

%% ----------------------------------------
%% TICKET
%% [prime_in_1 other_in_1] = scr_dsgn_preproc(proj,proj.param.n_trs_id1,run1_in_stim_times);
%% [prime_in_2 other_in_2] = scr_dsgn_preproc(proj,proj.param.n_trs_id2,run2_in_stim_times);
%% 
%% [prime_feel_1 other_feel_1] = scr_dsgn_preproc(proj,proj.param.n_trs_id1,run1_feel_stim_times);
%% [prime_feel_2 other_feel_2] = scr_dsgn_preproc(proj,proj.param.n_trs_id2,run2_feel_stim_times);

%% ----------------------------------------
%% TICKET
%% Should rewrite this whole system in terms of tasks and runs
%% to make it more flexible and general

%% load subjs
subjs = proj.process.subjs; 

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
    ex_betas = struct();

    %% -----------------------------------------
    %% LSA of scr signal (not preferred to LSS)
    %% path = [proj.path.physio.scr_clean,subj_study,'_',name,'_Identify_run_1.mat'];
    %% load(path);
    %% mdl_ex_1 = regstats(scr,prime_ex_1');
    %% betas1 = zscore(mdl_ex_1.beta(2:end)');
    %% 
    %% path = [proj.path.physio.scr_clean,subj_study,'_',name,'_Identify_run_2.mat'];
    %% load(path);
    %% mdl_ex_2 = regstats(scr,prime_ex_2');
    %% betas2 = zscore(mdl_ex_1.beta(2:end)');
    %% 
    %% grp_betas = [grp_betas;[betas1,betas2]];

    %% ----------------------------------------
    %% LSS of scr signal (Mumford, 2012) - Identify 1
    ex_betas.id1 = [];
    try
        path = [proj.path.physio.scr_clean,subj_study,'_',name,'_Identify_run_1.mat'];
        load(path);

        for j=1:size(prime_ex_1)
            prime = prime_ex_1(j,:);
            other = other_ex_1(j,:);
            mdl_ex_1 = regstats(scr,[prime_ex_1(j,:)',other_ex_1(j,:)']);
            ex_betas.id1 = [ex_betas.id1,mdl_ex_1.beta(2)'];
        end

        %%Normalize
        ex_betas.id1 = zscore(ex_betas.id1);

    catch
        logger(['  -LSS Error: SCR of Identify run 1: ',path],proj.path.logfile);
    end

    %% ----------------------------------------
    %% LSS of scr signal (Mumford, 2012) - Identify 2
    ex_betas.id2 = [];
    try
        path = [proj.path.physio.scr_clean,subj_study,'_',name,'_Identify_run_2.mat'];
        load(path);

        for j=1:size(prime_ex_2)
            prime = prime_ex_2(j,:);
            other = other_ex_2(j,:);
            mdl_ex_2 = regstats(scr,[prime_ex_2(j,:)',other_ex_2(j,:)']);
            ex_betas.id2 = [ex_betas.id2,mdl_ex_2.beta(2)'];
        end
        
       %%Normalize
       ex_betas.id2 = zscore(ex_betas.id2);

    catch
        logger(['  -LSS Error: SCR of Identify run 2: ',path],proj.path.logfile);
    end

    %% ----------------------------------------
    %% SAVE Individual Betas
    save([proj.path.betas.scr_ex_beta,subj_study,'_',name,'_ex_betas.mat'],'ex_betas');

    % debug
    if(~isempty(ex_betas.id1) & ~isempty(ex_betas.id2))
        ex_betas = [ex_betas.id1,ex_betas.id2];
        grp_betas = [grp_betas;ex_betas];
    end

end

%%Indicate completion of this process
proj.process.beta_scr_ex_id = 1;

%% Write out amended project params
save('proj.mat');

%% % debug
%% load(['/home/kabush/workspace/data/CTM/analysis/univ_lss_trgs/stim_a_scores.txt']);
%% load(['/home/kabush/workspace/data/CTM/analysis/univ_lss_trgs/' ...
%%       'stim_ids.txt']);
%% ex_ids = find(stim_ids==1);
%% grp_b = [];
%% for i=1:size(grp_betas,1)
%%     [b stat] = robustfit(grp_betas(i,:),stim_a_scores(ex_ids));
%%     grp_b = [grp_b,b(2)];
%% end
%% 
%% [h p ci stat] = ttest(grp_b);
%% p 
%% ci
