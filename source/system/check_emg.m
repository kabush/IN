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

logger('----------------------------------------',proj.path.logfile);
logger(' EMG Quality Check ',proj.path.logfile);
logger('----------------------------------------',proj.path.logfile);

%% Create the subjects to be analyzed (possible multiple studies)
subjs = proj.process.subjs; 

%% Preprocess fMRI of each subject in subjects list 
for i=1:numel(subjs)
    
    subj = subjs{i};

    %% extract subject info
    subj_study = subj.study;
    name = subj.name;

    % preprocessed files exist?
    subj.emg.ok = 0;

    subj.emg.id1.ok = 0;
    subj.emg.id2.ok = 0;

    subj.emg.corr.id1.ok = 0;
    subj.emg.corr.id2.ok = 0;
    subj.emg.corr.id1.exist = 0;
    subj.emg.corr.id1.non_zero = 0;
    subj.emg.corr.id2.exist = 0;
    subj.emg.corr.id2.non_zero = 0;

    subj.emg.zygo.id1.ok = 0;
    subj.emg.zygo.id2.ok = 0;
    subj.emg.zygo.id1.exist = 0;
    subj.emg.zygo.id1.non_zero = 0;
    subj.emg.zygo.id2.exist = 0;
    subj.emg.zygo.id2.non_zero = 0;
    
    %% Check for EMG corrugator files
    %  (identify 1)
    path = [proj.path.physio.emg_clean,subj_study,'_',name,'_Identify_run_1_corr.mat'];
    if exist(path,'file')==2
        subj.emg.corr.id1.exist = 1;
        if(subj.emg.corr.id1.exist)
            load(path);
            if(mean(rect_corr)>0)
                subj.emg.corr.id1.non_zero=1;
            end
        end
    end
    
    %  (identify 2)
    path = [proj.path.physio.emg_clean,subj_study,'_',name,'_Identify_run_2_corr.mat'];
    if exist(path,'file')==2
        subj.emg.corr.id2.exist = 1;
        if(subj.emg.corr.id2.exist)
            load(path);
            if(mean(rect_corr)>0)
                subj.emg.corr.id2.non_zero=1;
            end
        end
    end
    

    %% Check for EMG zygomaticus files
    %  (identify 1)
    path = [proj.path.physio.emg_clean,subj_study,'_',name,'_Identify_run_1_zygo.mat'];
    if exist(path,'file')==2
        subj.emg.zygo.id1.exist = 1;
        if(subj.emg.zygo.id1.exist)
            load(path);
            if(mean(rect_zygo)>0)
                subj.emg.zygo.id1.non_zero=1;
            end
        end
    end
    
    %  (identify 2)
    path = [proj.path.physio.emg_clean,subj_study,'_',name,'_Identify_run_2_zygo.mat'];
    if exist(path,'file')==2
        subj.emg.zygo.id2.exist = 1;
        if(subj.emg.zygo.id2.exist)
            load(path);
            if(mean(rect_zygo)>0)
                subj.emg.zygo.id2.non_zero=1;
            end
        end
    end
    
    %% ----------------------------------------
    %% Processing logic
    
    % check that all EMG corrugator  components are available (2 pieces)
    if(subj.emg.corr.id1.exist & subj.emg.corr.id1.non_zero)
        subj.emg.corr.id1.ok = 1;
    end

    if(subj.emg.corr.id2.exist & subj.emg.corr.id2.non_zero)
        subj.emg.corr.id2.ok = 1;
    end

    % check that all EMG zygomaticus  components are available (2 pieces)
    if(subj.emg.zygo.id1.exist & subj.emg.zygo.id1.non_zero)
        subj.emg.zygo.id1.ok = 1;
    end

    if(subj.emg.zygo.id2.exist & subj.emg.zygo.id2.non_zero)
        subj.emg.zygo.id2.ok = 1;
    end

    % combine recording types
    if(subj.emg.corr.id1.ok & subj.emg.zygo.id1.ok)
        subj.emg.id1.ok = 1;
    end

    if(subj.emg.corr.id2.ok & subj.emg.zygo.id2.ok)
        subj.emg.id2.ok = 1;
    end

    % combine runs
    if(subj.emg.id1.ok & subj.emg.id2.ok)
        subj.emg.ok = 1;
    else
        disp([subj_study,'_',name]);
    end

    % output
    proj.process.subjs{i} = subj;
    
end

%% Indicate quality check has been completed
proj.check.emg = 1;

%% ----------------------------------------
%% Write out amended project structure
save('proj.mat','proj');
