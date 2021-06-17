%%========================================
%%========================================
%%
%% Keith Bush, PhD (2017)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

%% Load in path data
load('proj.mat');

%% Initialize log section
logger(['********************************************'],proj.path.logfile);
logger([' Analyzing GS Classification (Refit)        '],proj.path.logfile);
logger(['********************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.analysis.gs_cls_refit_perm]);
    eval(['! rm -rf ',proj.path.analysis.gs_cls_refit_perm]);
    disp(['Creating ',proj.path.analysis.gs_cls_refit_perm]);
    eval(['! mkdir ',proj.path.analysis.gs_cls_refit_perm]);
end

%% ----------------------------------------
%% Load label for number of targets
label_id = load([proj.path.trg.ex,'stim_ids.txt']);
ex_id = find(label_id == proj.param.trg.ex_id);
N_ex = numel(ex_id);

%% ----------------------------------------
%% load subjs
subjs = proj.process.subjs;

all_raw_acc_v = [];
all_raw_acc_a = [];

all_raw_prm_v = [];
all_raw_prm_a = [];

good_sbj_ids = [];
for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;
    
    disp([subj_study,'_',name]);

    try

        % Load permutations
        load([proj.path.mvpa.fmri_ex_gs_perm,subj_study,'_',name,'_prds.mat']);
        prms = prds;

        % Load classification performance
        load([proj.path.mvpa.fmri_ex_gs_cls,subj_study,'_',name,'_prds.mat']);
        
        % Identify usable ids
        good_ids = ex_id; 

        % ----------------------------------------
        % GROUP-LEVEL PERMUTATION STAT
        %
        % Extract mean accuracy for this subject on usable ids
        acc_v = zeros(1,N_ex)-1;
        acc_a = zeros(1,N_ex)-1;

        perm_v = zeros(1,N_ex)-1;
        perm_a = zeros(1,N_ex)-1;

        for j = 1:numel(good_ids)

            id = good_ids(j); 
            indx = find(ex_id==id);

            acc_v(indx) = mean(prds.v_cls_acc(:,j),1);
            acc_a(indx) = mean(prds.a_cls_acc(:,j),1);

            %% Group permutation values (for group comparison)
            prm_v(indx) = mean(prms.v_cls_perm(:,j),1);
            prm_a(indx) = mean(prms.a_cls_perm(:,j),1);

        end

        % Take column means (explicit 1 to handle nrows=1 case)
        all_raw_acc_v = [all_raw_acc_v;acc_v];
        all_raw_acc_a = [all_raw_acc_a;acc_a];

        all_raw_prm_v = [all_raw_prm_v;prm_v];
        all_raw_prm_a = [all_raw_prm_a;prm_a];

        % ----------------------------------------
        % SUBJECT-LEVEL PERMUTATION STAT
        %
        % Compile permutation test for single subject test
        smp_prm_v = [];
        smp_prm_a = [];
        
        Nperm = 1000;
        Nresample = 30;
        for k = 1:Nperm
        
            perm_smp_v = [];
            perm_smp_a = [];
        
            perm_ids = randperm(size(prms.v_cls_perm,1));
            perm_ids = perm_ids(1:Nresample);
        
            % Construct the row
            for j = 1:numel(good_ids)
                perm_smp_v = [perm_smp_v,mean(prms.v_cls_perm(perm_ids,j),1)];
                perm_smp_a = [perm_smp_a,mean(prms.a_cls_perm(perm_ids,j),1)];
            end
        
            smp_prm_v = [smp_prm_v;perm_smp_v];
            smp_prm_a = [smp_prm_a;perm_smp_a];
        
        end

       
        sbj_prms{i}.smp_prm_v = smp_prm_v;
        sbj_prms{i}.smp_prm_a = smp_prm_a;

        good_sbj_ids = [good_sbj_ids,i];
        
    catch
        disp(' no beta series');
    end

end


%% ----------------------------------------
%% "Knows-what-it-knows" performance

%% VALENCE Analysis
logger(['****VALENCE Analysis****'],proj.path.logfile);
calc_mvpa_ex_gs_refit_perm(proj,good_sbj_ids,all_raw_acc_v,all_raw_prm_v,sbj_prms,'v');

%% AROUSAL Analysis
logger(['****AROUSAL Analysis****'],proj.path.logfile);
calc_mvpa_ex_gs_refit_perm(proj,good_sbj_ids,all_raw_acc_a,all_raw_prm_a,sbj_prms,'a');


