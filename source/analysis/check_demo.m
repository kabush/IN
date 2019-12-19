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

%% Initialize log section
logger(['************************************************'],proj.path.logfile);
logger([' Characterizing Age|Sex of Subjects             '],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

m_age = [];
f_age = [];


%% ----------------------------------------
%% iterate over study subjects
for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    %% debug
    logger([subj_study,':',name],proj.path.logfile);

    try

        %% Load subject's sex
        demo = readtable([proj.path.raw_data,proj.path.demo,'/',subj_study,'.csv']);
        id = find(strcmp(demo.ID,name)~=0);
        sex = demo.Type(id);
        age = demo.Age(id)

        if(sex==1)
            m_age = [m_age,age];
        else
            f_age = [f_age,age];
        end

    catch 
        disp(['  Error: could not load file(s)']);
    end

end

