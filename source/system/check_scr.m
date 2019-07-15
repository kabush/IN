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
logger(' SCR Quality Check ');
logger('----------------------------------------');

%% Create the subjects to be analyzed (possible multiple studies)
subjs = proj.process.subjs;

%% Preprocess fMRI of each subject in subjects list 
for i=1:numel(subjs)
    
    subj = subjs{i};

    %% extract subject info
    subj_study = subj.study;
    name = subj.name;

    % preprocessed files exist?
    subj.scr.ok = 0;
    subj.scr.id1.ok = 0;
    subj.scr.id2.ok = 0;
    
    % check for existence of clean scr file (identify 1)
    if exist([proj.path.physio.scr_clean,subj_study,'_',name,'_Identify_run_1.mat'],'file')==2
        subj.scr.id1.ok = 1;
    end
    
    % check for existence of clean scr file (identify 1)
    if exist([proj.path.physio.scr_clean,subj_study,'_',name,'_Identify_run_2.mat'],'file')==2
        subj.scr.id2.ok = 1;
    end
    
    %% ----------------------------------------
    %% Processing logic
    
    % check that all fmri components are available (4 pieces)
    if(subj.scr.id1.ok+subj.scr.id2.ok==2);
        subj.scr.ok = 1;
    else
        disp([subj_study,'_',name]);
    end
    
    % *** TICEKT ***don't exclude based on physio (secondary measure, handle
    % in specific analyses
    % % assign master flag to false 
    % if(subj.scr.ok==0)
    %     subj.ok = 0;
    % end
    
    % output
    proj.process.subjs{i} = subj;
    
end

%% Indicate quality check has been completed
proj.check.scr = 1;

%% ----------------------------------------
%% Write out amended project structure
save('proj.mat','proj');
