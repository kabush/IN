%%========================================
%%========================================
%%
%% Keith Bush, PhD (2019)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

%% Load in path data
load('proj.mat');

%% Set-up Directory Structure
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.physio.emg_clean]);
    eval(['! rm -rf ',proj.path.physio.emg_clean]);
    disp(['Creating ',proj.path.physio.emg_clean]);
    eval(['! mkdir ',proj.path.physio.emg_clean]);
end

%% Create the subjects to be analyzed (possible multiple studies)
subjs = load_subjs(proj);
disp(['Processing EMG of ',num2str(numel(subjs)),' subjects']);

%% ----------------------------------------
%% Compute EMGs over each subject individually
for i=1:numel(subjs)
   
    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    %% debug
    disp([subj_study,':',name]);

    try

        %% ----------------------------------------
        %% Process Identify Run 1
        disp('Processing (in Matlab) Identify 1');
        
        %% load Identify 1 data
        in_path = [proj.path.raw_data,subj_study,'/', ...
                   proj.path.raw_physio,'/',subj_study,'_',name,'/', ...
                   subj_study,'_',name,'_Identify_run_1.mat'];
        load(in_path);
   
        size(data)
     
        % subselect component from BIOPAC
        chan_zygo = proj.param.physio.chan_emg_zygo;
        chan_corr = proj.param.physio.chan_emg_corr;
        Ntrs = proj.param.mri.n_trs_id1;    
        data_zygo = data(:,chan_zygo);
        data_corr = data(:,chan_corr);
        
        %process zygo
        emg_zygo = emg_preproc(proj,Ntrs,data_zygo);
        rect_zygo = abs(emg_zygo);
        
        %process corr
        emg_corr = emg_preproc(proj,Ntrs,data_corr);
        rect_corr = abs(emg_corr);
        
        %save emg
        out_path = [proj.path.physio.emg_clean,subj_study,'_',name];
        save([out_path,'_run_1_corr.mat'],'rect_corr');
        save([out_path,'_run_1_zygo.mat'],'rect_zygo');
        
        %% ----------------------------------------
        %% Process Identify Run 2
        disp('Processing (in Matlab) Identify 2');
        
        %% load Identify 1 data
        in_path = [proj.path.raw_data,subj_study,'/', ...
                   proj.path.raw_physio,'/',subj_study,'_',name,'/', ...
                   subj_study,'_',name,'_Identify_run_2.mat'];
        load(in_path);
        
        % subselect component from BIOPAC
        chan_zygo = proj.param.physio.chan_emg_zygo;
        chan_corr = proj.param.physio.chan_emg_corr;
        Ntrs = proj.param.mri.n_trs_id2;    
        data_zygo = data(:,chan_zygo);
        data_corr = data(:,chan_corr);
        
        %process zygo
        emg_zygo = emg_preproc(proj,Ntrs,data_zygo);
        rect_zygo = abs(emg_zygo);
        
        %process corr
        emg_corr = emg_preproc(proj,Ntrs,data_corr);
        rect_corr = abs(emg_corr);
        
        %save emg
        out_path = [proj.path.physio.emg_clean,subj_study,'_',name];
        save([out_path,'_run_2_corr.mat'],'rect_corr');
        save([out_path,'_run_2_zygo.mat'],'rect_zygo');
        
    catch
        disp('   EMG data may be unavailable');
    end
    
end
