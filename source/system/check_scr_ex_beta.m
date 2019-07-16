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
logger(' SCR Beta-series Quality Check ');
logger('----------------------------------------');

%% Create the subjects to be analyzed (possible multiple studies)
subjs = proj.process.subjs; %load_subjs(proj);

%% Preprocess fMRI of each subject in subjects list 
for i=1:numel(subjs)
    
    subj = subjs{i};
    
    %% extract subject info
    subj_study = subj.study;
    name = subj.name;
    
    % files are usable quality
    subj.beta.scr_ex_id.ok = 0; %overall usable
    
    subj.beta.scr_ex_id1.ok = 0; %usable
    subj.beta.scr_ex_id1.exist = 0; %exist
    subj.beta.scr_ex_id1.nan_ok = 0; %exist
    subj.beta.scr_ex_id1.nan_ids = [];
    subj.beta.scr_ex_id2.ok = 0; %usable
    subj.beta.scr_ex_id2.exist = 0; %exist
    subj.beta.scr_ex_id2.nan_ok = 0; %exist
    subj.beta.scr_ex_id2.nan_ids = [];
    
    % check for existence of Identify 1 beta-series
    path = [proj.path.betas.scr_ex_beta,subj_study,'_',name,'_ex_betas.mat'];
    if exist(path,'file')==2
        load(path);
        
        %% check existence: id 1
        if(numel(ex_betas.id1)>0)
            subj.beta.scr_ex_id1.exist = 1;
        else
            disp([subj_study,'_',name,': ID 1 ~exist']);
        end
        
        %% check values
        if(subj.beta.scr_ex_id1.exist)
            if(numel(find(isnan(ex_betas.id1)))==0)
                subj.beta.scr_ex_id1.nan_ok = 1;
            else
                disp([subj_study,'_',name,': ID 1 nan']);
                subj.beta.scr_ex_id1.nan_ids = find(isnan(ex_betas.id1));
            end
        end
        
        %% check existence: id 2
        if(numel(ex_betas.id2)>0)
            subj.beta.scr_ex_id2.exist = 1;
        else
            disp([subj_study,'_',name,': ID 2 ~exist']);
        end
        
        %% check values
        if(subj.beta.scr_ex_id2.exist)
            if(numel(find(isnan(ex_betas.id2)))==0)
                subj.beta.scr_ex_id2.nan_ok = 1;
            else
                disp([subj_study,'_',name,': ID 2 nan']);
                subj.beta.scr_ex_id2.nan_ids = find(isnan(ex_betas.id2));
            end
        end
        
    end
    
    %% ----------------------------------------
    %% Quality control logic
    
    %% Check individual id betas
    if(subj.beta.scr_ex_id1.exist & subj.beta.scr_ex_id1.nan_ok)
        subj.beta.scr_ex_id1.ok = 1;
    end
    
    %% Check individual id betas
    if(subj.beta.scr_ex_id2.exist & subj.beta.scr_ex_id2.nan_ok)
        subj.beta.scr_ex_id2.ok = 1;
    end
    
    %% combine runs
    if(subj.beta.scr_ex_id1.ok & subj.beta.scr_ex_id2.ok)
        subj.beta.scr_ex_id.ok = 1;
    end
    
    % assign master flag to false
    if(subj.beta.scr_ex_id.ok==0)
        subj.ok = 0;
    end
    
    % output
    proj.process.subjs{i} = subj;
    
end

%% Indicate quality check has been completed
proj.check.beta_scr_ex_id = 1;

%% ----------------------------------------
%% Write out amended project structure
save('proj.mat','proj');