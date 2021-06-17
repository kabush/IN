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
logger(['Permutation of Intra-subject LOOCV GS MVPA '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.mvpa.fmri_ex_gs_perm]);
    eval(['! rm -rf ',proj.path.mvpa.fmri_ex_gs_perm]);
    disp(['Creating ',proj.path.mvpa.fmri_ex_gs_perm]);
    eval(['! mkdir ',proj.path.mvpa.fmri_ex_gs_perm]);
end

%% ----------------------------------------
%% Load labels;
v_label = load([proj.path.trg.ex,'stim_v_labs.txt']);
a_label = load([proj.path.trg.ex,'stim_a_labs.txt']);
label_id = load([proj.path.trg.ex,'stim_ids.txt']);
v_score = load([proj.path.trg.ex,'stim_v_scores.txt']);
a_score = load([proj.path.trg.ex,'stim_a_scores.txt']);

%% ----------------------------------------
%% load subjs
subjs = proj.process.subjs; 

%% ----------------------------------------
%% allocate storage
all_v_cls_acc = [];
all_a_cls_acc = [];

%% ----------------------------------------
%% iterate over study subjects
for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    %% log
    logger([subj_study,':',name],proj.path.logfile);

    try

        %% Load gray matter mask 
        gm_nii = load_nii([proj.path.mri.gm_mask,subj_study,'.',name,'.gm.nii']);
        mask = double(gm_nii.img);
        brain_size=size(mask);
        mask = reshape(mask,brain_size(1)*brain_size(2)*brain_size(3),1);
        in_brain=find(mask==1);  
        
        %% Load beta-series
        base_nii = load_nii([proj.path.betas.fmri_ex_beta,subj_study,'_',name,'_lss.nii']);
        brain_size = size(base_nii.img);
        
        %% Vectorize the base image
        base_img = vec_img_2d_nii(base_nii);
        base_img = reshape(base_img,brain_size(1)*brain_size(2)*brain_size(3),brain_size(4));
        
        %% Concatenate the MASKED base image
        subj_img = base_img(in_brain,:)';
        
        %% Concatenate all label/subj identifiers
        subj_id = repmat(id,numel(v_label),1);
        subj_i = repmat(i,numel(v_label),1);
        
        %% Subselect extrinsic data
        ex_id = find(label_id==proj.param.trg.ex_id);
        ex_id_good = ex_id; % over-ride this check 

        %% Log which ex_ids are used for model fitting
        proj.process.subjs{i}.beta.mri_ex_id.good_id = ex_id_good;

        ex_img = subj_img(ex_id_good,:);
        ex_subj_id = subj_id(ex_id_good,1);
        ex_v_label = v_label(ex_id_good,1);
        ex_a_label = a_label(ex_id_good,1);
        
        %% Perform Gram-Schmidt on data
        [~,ex_img] = gram_schmidt(ex_img);
        
        %% Peform quality check of generated features
        qlty = check_gs_img_qlty(ex_img);
        
        if(qlty.ok)
            
            %% Initialize the prediction structure of this subject
            prds = struct();
            prds.v_cls_perm = [];
            prds.a_cls_perm = [];

            %% ----------------------------------------
            %% Classify extrinsic VALENCE examples
            for j=1:proj.param.mvpa.n_perm

                rnd_ids = randperm(numel(ex_v_label));

                %% Fit the data space
                [~,~,v_tst_hd,~,v_cls_stats] = classify_loocv(ex_img,ex_v_label(rnd_ids),...
                                                              ex_subj_id,id, proj.param.mvpa.kernel);
                
                %% Store results
                v_cls_acc = cell2mat(v_cls_stats.tst_acc);
                prds.v_cls_perm = [prds.v_cls_perm;v_cls_acc];

            end
            
            % debug
            all_v_cls_acc = [all_v_cls_acc;mean(mean(prds.v_cls_perm,1))];
            logger(['  v acc: ',num2str(mean(mean(prds.v_cls_perm,2)))],proj.path.logfile);
            
            %% ----------------------------------------
            %% Classify extrinsic AROUSAL examples
            for j=1:proj.param.mvpa.n_perm
            
                rnd_ids = randperm(numel(ex_a_label));
            
                %% Fit the data space
                [~,~,a_tst_hd,~,a_cls_stats] = classify_loocv(ex_img,ex_a_label(rnd_ids),...
                                                              ex_subj_id,id, proj.param.mvpa.kernel);
                
                %% Store results
                a_cls_acc = cell2mat(a_cls_stats.tst_acc);
                prds.a_cls_perm = [prds.a_cls_perm;a_cls_acc];
            
            end
            
            % debug
            all_a_cls_acc = [all_a_cls_acc;mean(mean(prds.a_cls_perm,1))];
            logger(['  a acc: ',num2str(mean(mean(prds.a_cls_perm,2)))],proj.path.logfile);

            %% ----------------------------------------
            %% Save out results
            save([proj.path.mvpa.fmri_ex_gs_perm,subj_study,'_',name,'_prds.mat'],'prds');

        else
            logger(['  -MVPA Error: ex_img failed quality check'],proj.path.logfile);
        end

    catch
        logger(['  -MVPA Error: possible missing beta series'],proj.path.logfile);
    end

end

% log summary results
[h p ci stat] = ttest(all_v_cls_acc);
logger(['grp v acc ci=[',num2str(ci(1)),',',num2str(ci(2)),['], ' ...
                    'p='],num2str(p)],proj.path.logfile);

[h p ci stat] = ttest(all_a_cls_acc);
logger(['grp a acc ci=[',num2str(ci(1)),',',num2str(ci(2)),['], ' ...
                    'p='],num2str(p)],proj.path.logfile);
