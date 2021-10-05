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
logger(' Summarize Project Data',proj.path.logfile);
logger('----------------------------------------',proj.path.logfile);

%% Create the subjects to be analyzed (possible multiple studies)
subjs = proj.process.subjs; 

%% Preprocess analysis
mri_cnt = 0;
scr_cnt = 0;
emg_cnt = 0;

%% Beta-series analysis
beta_mri_ex_id_cnt = 0;
beta_scr_ex_id_cnt = 0;

%% Overall anlaysis
tot_cnt = 0;

for i=1:numel(subjs)

    subj = proj.process.subjs{i};

    if(proj.process.mri & proj.check.mri)
        mri_cnt = mri_cnt + subj.mri.ok;
    end

    if(proj.process.scr & proj.check.scr)
        scr_cnt = scr_cnt + subj.scr.ok;
    end

    if(proj.process.emg & proj.check.emg)
        emg_cnt = emg_cnt + subj.emg.ok;
    end

    if(proj.process.beta_mri_ex_id & proj.check.beta_mri_ex_id)
        beta_mri_ex_id_cnt = beta_mri_ex_id_cnt + subj.beta.mri_ex_id.ok;
    end

    if(proj.process.beta_scr_ex_id & proj.check.beta_scr_ex_id)
        beta_scr_ex_id_cnt = beta_scr_ex_id_cnt + subj.beta.scr_ex_id.ok;
    end

    if(proj.process.mri & proj.check.mri & proj.process.scr & ...
       proj.process.emg & proj.process.beta_mri_ex_id & ...
       proj.process.beta_scr_ex_id & proj.check.mri & proj.check.scr ...
       & proj.check.emg & proj.check.beta_mri_ex_id & ...
       proj.check.beta_scr_ex_id)
        if(subj.mri.ok+subj.scr.ok+subj.emg.ok+ ...
           subj.beta.mri_ex_id.ok+subj.beta.scr_ex_id.ok==5)
            tot_cnt = tot_cnt + 1;
        else
            disp(['Missing Data: ',subj.study,'_',subj.name]);
        end
    end

end

logger(['--Passed OVERALL qlty   ',num2str(100*(tot_cnt/numel(subjs))),'%'],proj.path.logfile);
logger(['--Passed MRI qlty:      ',num2str(100*(mri_cnt/numel(subjs))),'%'],proj.path.logfile);
logger(['--Passed SCR qlty:      ',num2str(100*(scr_cnt/numel(subjs))),'%'],proj.path.logfile);
logger(['--Passed EMG qlty:      ',num2str(100*(emg_cnt/numel(subjs))),'%'],proj.path.logfile);
logger(['--Passed MRI beta qlty: ',num2str(100*(beta_mri_ex_id_cnt/numel(subjs))),'%'],proj.path.logfile);
logger(['--Passed SCR beta qlty: ',num2str(100*(beta_scr_ex_id_cnt/numel(subjs))),'%'],proj.path.logfile);