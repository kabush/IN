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

%% Set-up Directory Structure
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.mri.mri_clean]);
    eval(['! rm -rf ',proj.path.mri.mri_clean]);
    disp(['Creating ',proj.path.mri.mri_clean]);
    eval(['! mkdir ',proj.path.mri.mri_clean]);
end

%% Create the subjects to be analyzed (possible multiple studies)
subjs = load_subjs(proj);

logger(['************************************'],proj.path.logfile);
logger(['Prepocessing fMRI of ',num2str(numel(subjs)),' subjects'],proj.path.logfile);
logger(['************************************'],proj.path.logfile);

%% Preprocess fMRI of each subject in subjects list 
for i=1:numel(subjs)

    %% Get data from raw and formatted into the correct name convention
    dlmwrite([proj.path.code,'tmp/home_path.txt'],proj.path.home,'');
    dlmwrite([proj.path.code,'tmp/raw_path.txt'],proj.path.raw_data,'');
    dlmwrite([proj.path.code,'tmp/lib_path.txt'],proj.path.tools.kablab,'');
    dlmwrite([proj.path.code,'tmp/project_name.txt'],proj.path.name,'');

    %% Get preprocessing params
    dlmwrite([proj.path.code,'tmp/tr.txt'],num2str(proj.param.mri.TR),'');
    dlmwrite([proj.path.code,'tmp/slices.txt'],proj.param.mri.slices,'');
    dlmwrite([proj.path.code,'tmp/slice_pattern.txt'],proj.param.mri.slice_pattern,'');
    dlmwrite([proj.path.code,'tmp/do_anat.txt'],proj.param.mri.do_anat,'');
    dlmwrite([proj.path.code,'tmp/do_epi.txt'],proj.param.mri.do_epi,'');
    dlmwrite([proj.path.code,'tmp/tasks.txt'],proj.param.mri.tasks,'');
    dlmwrite([proj.path.code,'tmp/scans.txt'],proj.param.mri.scans,'');
    dlmwrite([proj.path.code,'tmp/rest_scans.txt'],proj.param.mri.rest_scans,'');

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    %% debug
    logger([subj_study,':',name],proj.path.logfile);

    %% Get data from raw and formatted into the correct name convention
    dlmwrite([proj.path.code,'tmp/study.txt'],subj_study,'');
    dlmwrite([proj.path.code,'tmp/subject.txt'],name,'');

    %% Do the preprocessing
    eval(['! ',proj.path.code,'source/preprocess/preprocess_fmri_afni ',proj.path.home,' ',proj.path.name]);

    %% Clean-up
    eval(['! rm ',proj.path.code,'tmp/*']);

end

%%Indicate completion of this process
proj.process.mri = 1;

%% Write out amended project params
save('proj.mat');
