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
logger(['*************************************************'],proj.path.logfile);
logger(['Analyzing ER Skill (VALENCE SEX-DIFFS)           '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.analysis.vr_skill]);
    eval(['! rm -rf ',proj.path.analysis.vr_skill]);
    disp(['Creating ',proj.path.analysis.vr_skill]);
    eval(['! mkdir ',proj.path.analysis.vr_skill]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% ----------------------------------------
%% load subjs

b_all = []; %% for power
sex = [];

for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    % log analysis of subject
    logger([subj_study,'_',name],proj.path.logfile);

    try

        %% Load IN trajectory structures
        load([proj.path.ctrl.in_dyn,subj_study,'_',name, ...
              '_prds.mat']);

        if(isfield(prds,'v_dcmp'))
            
            %% extract stims and mean "feel"
            raw_stim = prds.v_dcmp.stim;
            [~,ids] = sort(raw_stim');
            stim = raw_stim(ids,:);
            feel = mean(prds.v_dcmp.feel(ids,:),2);

            %% remove extreme outliers (TICKET hardcoded outliers)
            stim_keep_ids = find(abs(stim)<=3);
            stim_feel_ids = find(abs(feel)<=3);
            cmb_keep_ids = intersect(stim_keep_ids,stim_feel_ids);
            disp(['                  ', num2str(numel(cmb_keep_ids))]);
            stim_clean = stim(cmb_keep_ids);
            feel_clean = feel(cmb_keep_ids);
         
            [b stat] = robustfit(stim_clean,feel_clean);
            b_all = [b_all,b(2)]; %% for power

            %% ****************************************
            %% GET SEX OF SUBJ
            %% ****************************************
            demo = readtable(['/raw/bush/demo/',subj_study,'.csv']);
            id = find(strcmp(demo.ID,name)~=0);
            sex = [sex,demo.Type(id)];
            
        else
            disp(['  -Could not find v_dcmp for: ',subj_study,'_', ...
                  name],proj.path.logfile);
        end
        
    catch
        % do nothing
        logger(['  -Could not find/load prds for: ',subj_study,'_', ...
                name],proj.path.logfile);
    end

end


b_m = b_all(find(sex==1));
b_f = b_all(find(sex==2));

p = ranksum(b_m,b_f);
disp(['ranksum of b_all sex diffs, p=',num2str(p)]);
