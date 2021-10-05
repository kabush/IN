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
logger([' Global Permutation Testing of MVPA Hyperplanes '],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.haufe.fmri_ex_gm_mdl]);
    eval(['! rm -rf ',proj.path.haufe.fmri_ex_gm_mdl]);
    disp(['Creating ',proj.path.haufe.fmri_ex_gm_mdl]);
    eval(['! mkdir ',proj.path.haufe.fmri_ex_gm_mdl]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% ----------------------------------------
%% load labels
v_label = load([proj.path.trg.ex,'stim_v_labs.txt']);
a_label = load([proj.path.trg.ex,'stim_a_labs.txt']);
label_id = load([proj.path.trg.ex,'stim_ids.txt']);

%% ----------------------------------------
%% iterate over permuations
Nperm = proj.param.haufe.npermute;
Nloop = Nperm + 1; %% first loop is true model structure
Nchunk = proj.param.haufe.chunk;

%% storage for group Haufe
grp_haufe_v = zeros(172800,Nloop);
grp_haufe_a = zeros(172800,Nloop);

%% permutation significance levels
alpha05 = 0.05;
alpha01 = 0.01;
alpha001 = 0.001;

for i = 1:Nloop

    %%storage for group haufe 
    all_haufe_v_wts = zeros(172800,numel(subjs));
    all_haufe_v_mask = zeros(172800,numel(subjs));

    all_haufe_a_wts = zeros(172800,numel(subjs));
    all_haufe_a_mask = zeros(172800,numel(subjs));
    
    qlty_vec = zeros(1,numel(subjs));

    %% ----------------------------------------
    %% iterate over study subjects
    for j = 1:numel(subjs)
        
        %% extract subject info
        subj_study = subjs{j}.study;
        name = subjs{j}.name;
        id = subjs{j}.id;
        
        %% debug
        logger([subj_study,':',name,':i=',num2str(i)],proj.path.logfile);
        
        try
            
            %% Load gray matter mask 
            gm_nii = load_untouch_nii([proj.path.mri.gm_mask,subj_study,'.',name,'.gm.nii']);
            mask = double(gm_nii.img);
            brain_size=size(mask);
            mask = reshape(mask,brain_size(1)*brain_size(2)*brain_size(3),1);
            in_brain=find(mask==1);  
            
            %% Load beta-series
            base_nii = load_untouch_nii([proj.path.betas.fmri_ex_beta,subj_study,'_',name,'_lss.nii']);
            brain_size = size(base_nii.img);
            
            %% Vectorize the base image
            base_img = vec_img_2d_nii(base_nii);
            base_img = reshape(base_img,brain_size(1)*brain_size(2)*brain_size(3),brain_size(4));
            
            %% Concatenate the MASKED base image
            all_img = base_img(in_brain,:)';
            
            %% Concatenate all label/subj identifiers
            subj_id = repmat(id,numel(v_label),1);
            subj_j = repmat(j,numel(v_label),1);
            
            %% Subselect extrinsic data
            ex_id = find(label_id==proj.param.trg.ex_id);
            ex_img = all_img(ex_id,:);
            ex_subj_id = subj_id(ex_id,1);
            ex_v_label = v_label(ex_id,1);
            ex_a_label = a_label(ex_id,1);
            
            %% Peform quality check of generated features
            qlty = check_gm_img_qlty(ex_img);
            
            if(qlty.ok)
                
                %% ----------------------------------------
                %% VALENCE MODEL

                %% Balance positive and negative VALENCE examples
                v_pos_ids = find(ex_v_label==1);
                v_neg_ids = find(ex_v_label==-1);
                Npos = numel(v_pos_ids);
                Nneg = numel(v_neg_ids);
                
                if(Npos >= Nneg)
                    Nsample = Nneg;
                else
                    Nsample = Npos;
                end
                
                %% Randomly re-order and combine samples
                rnd_v_pos_ids = v_pos_ids(randsample(1:numel(v_pos_ids),Nsample))';
                rnd_v_neg_ids = v_neg_ids(randsample(1:numel(v_neg_ids),Nsample))';
                rnd_v_cmb_ids = [rnd_v_pos_ids,rnd_v_neg_ids];
                
                %% Reassign ids separately for features and labels
                feat_v_ids = rnd_v_cmb_ids;
                lab_v_ids = rnd_v_cmb_ids; 

                %% Only first iteration is structure (remaining
                %% loops are permutations, therefore randomize labels
                if(i>1)
                    lab_v_ids = randsample(lab_v_ids,numel(lab_v_ids));
                end

                %% Fit classifier
                v_model = fitcsvm(ex_img(feat_v_ids,:),ex_v_label(lab_v_ids,1), ...
                                  'KernelFunction',proj.param.mvpa.kernel);


                %% Construct Valence Haufe tranform
                wts = v_model.Beta;
                haufe_v_wts = zscore(fast_haufe(ex_img(feat_v_ids,:),wts,Nchunk));
                all_haufe_v_wts(in_brain,j) = haufe_v_wts;
                all_haufe_v_mask(in_brain,j) = 1;
                
                %% ----------------------------------------
                %% AROUSAL MODEL
                
                %% Balance positive and negative VALENCE examples
                a_pos_ids = find(ex_a_label==1);
                a_neg_ids = find(ex_a_label==-1);
                Npos = numel(a_pos_ids);
                Nneg = numel(a_neg_ids);
                
                if(Npos >= Nneg)
                    Nsample = Nneg;
                else
                    Nsample = Npos;
                end
                
                %% Randomly re-order and combined samples
                rnd_a_pos_ids = a_pos_ids(randsample(1:numel(a_pos_ids),Nsample))';
                rnd_a_neg_ids = a_neg_ids(randsample(1:numel(a_neg_ids),Nsample))';
                rnd_a_cmb_ids = [rnd_a_pos_ids,rnd_a_neg_ids];

                %% Reassign ids separately for features and labels
                feat_a_ids = rnd_a_cmb_ids;
                lab_a_ids = rnd_a_cmb_ids; 

                %% Only first iteration is structure (remaining
                %% loops are permutations, therefore randomize labels
                if(i>1)
                    lab_a_ids = randsample(lab_a_ids,numel(lab_a_ids));
                end

                %% Fit classifier
                a_model = fitcsvm(ex_img(feat_a_ids,:),ex_a_label(lab_a_ids,1), ...
                                  'KernelFunction',proj.param.mvpa.kernel);

                %% Construct Arousal Haufe transform
                wts = a_model.Beta;
                haufe_a_wts = zscore(fast_haufe(ex_img(feat_a_ids,:),wts,Nchunk));
                all_haufe_a_wts(in_brain,j) = haufe_a_wts;
                all_haufe_a_mask(in_brain,j) = 1;

                qlty_vec(j)=1;
                
            end
            
        catch
            logger(['  -Haufe Error'],proj.path.logfile);
        end

    end

    %% ----------------------------------------
    %% Extract Quality fits        
    qlty_ids = find(qlty_vec==1);
    qlty_n = numel(qlty_ids);
    
    qlty_haufe_v_wts = all_haufe_v_wts(:,qlty_ids);
    qlty_haufe_v_mask = all_haufe_v_mask(:,qlty_ids);
    
    qlty_haufe_a_wts = all_haufe_a_wts(:,qlty_ids);
    qlty_haufe_a_mask = all_haufe_a_mask(:,qlty_ids);
    
    %% ----------------------------------------
    %% Group Mean Haufe transforms (1 saved per permutation)
    
    %% Group valence Haufe
    ahf_sum = sum(qlty_haufe_v_mask,2);
    row_ids_v = find(ahf_sum>(qlty_n/2));
    grp_haufe_v(row_ids_v,i) = mean(qlty_haufe_v_wts(row_ids_v,:),2);

    %% T-score version
    if(i==1)
        grp_haufe_v_tstat = 0*grp_haufe_v(:,1);
        for k=1:numel(row_ids_v)
            [h p ci stat] = ttest(qlty_haufe_v_wts(row_ids_v(k),:));
            grp_haufe_v_tstat(row_ids_v(k),1)=stat.tstat;
        end
    end

    
    %% Group arousal Haufe
    ahf_sum = sum(qlty_haufe_a_mask,2);
    row_ids_a = find(ahf_sum>(qlty_n/2));
    grp_haufe_a(row_ids_a,i) = mean(qlty_haufe_a_wts(row_ids_a,:),2);
    
    %% T-score version
    if(i==1)
        grp_haufe_a_tstat = 0*grp_haufe_a(:,1);
        for k=1:numel(row_ids_v)
            [h p ci stat] = ttest(qlty_haufe_a_wts(row_ids_a(k),:));
            grp_haufe_a_tstat(row_ids_a(k),1)=stat.tstat;
        end
    end

    %% ----------------------------------------
    %% Do permutation test given samples available
    
    if(i>1)

        save([proj.path.haufe.fmri_ex_gm_mdl,'grp_haufe_v_n=',num2str(i-1),'_of_N=',num2str(Nperm),'.mat'],'grp_haufe_v');
        save([proj.path.haufe.fmri_ex_gm_mdl,'grp_haufe_a_n=',num2str(i-1),'_of_N=',num2str(Nperm),'.mat'],'grp_haufe_a');

        if(i>2)
            eval(['! rm ',proj.path.haufe.fmri_ex_gm_mdl,'grp_haufe_v_n=',num2str(i-2),'_of_N=',num2str(Nperm),'.mat']);
            eval(['! rm ',proj.path.haufe.fmri_ex_gm_mdl,'grp_haufe_a_n=',num2str(i-2),'_of_N=',num2str(Nperm),'.mat']);
        end
    
        %% ----------------------------------------
        %% Valence
        sig_ids_05_v = [];
        sig_ids_01_v = [];
        sig_ids_001_v = [];
        
        for k=1:numel(row_ids_v)
            
            % Count extrem random samples
            Next = 0;
            if(grp_haufe_v(row_ids_v(k),1)>0)
                Next = numel(find(grp_haufe_v(row_ids_v(k),2:i)>grp_haufe_v(row_ids_v(k),1)));
            else
                Next = numel(find(grp_haufe_v(row_ids_v(k),2:i)<grp_haufe_v(row_ids_v(k),1)));
            end
            
            % Do 2-sided tests
            if(Next<round((alpha05/2)*i))
                sig_ids_05_v = [sig_ids_05_v,row_ids_v(k)];
            end
            
            if(Next<round((alpha01/2)*i))
                sig_ids_01_v = [sig_ids_01_v,row_ids_v(k)];
            end
            
            if(Next<round((alpha001/2)*i))
                sig_ids_001_v = [sig_ids_001_v,row_ids_v(k)];
            end
            
        end

        % ----------------------------------------
        % Save out: mean encoding of group gray-matter voxels
        if(numel(row_ids_v)>0)

            % mu_v_haufe_nii = build_nii_from_gm_mask(grp_haufe_v(row_ids_v,1),gm_nii,row_ids_v);
            % save_nii(mu_v_haufe_nii,[proj.path.haufe.fmri_ex_gm_mdl,'mu_haufe_v_N=',num2str(Nperm),'.nii']);

            mu_v_haufe_nii = build_nii_from_gm_mask(grp_haufe_v_tstat(row_ids_v,1),gm_nii,row_ids_v);
            save_untouch_nii(mu_v_haufe_nii,[proj.path.haufe.fmri_ex_gm_mdl,'mu_haufe_v_N=',num2str(Nperm),'.nii']);

        end

        % ----------------------------------------
        % Save out: mean encoding of permstrap sign. (p<0.05) group
        % gray-matter voxels
        if(numel(sig_ids_05_v)>0)

            % mu_perm_v_haufe_nii = build_nii_from_gm_mask(grp_haufe_v(sig_ids_05_v,1),gm_nii,sig_ids_05_v);
            % save_nii(mu_perm_v_haufe_nii,[proj.path.haufe.fmri_ex_gm_mdl,'mu_perm_haufe_v_N=',num2str(Nperm),'_05.nii']);

            mu_perm_v_haufe_nii = build_nii_from_gm_mask(grp_haufe_v_tstat(sig_ids_05_v,1),gm_nii,sig_ids_05_v);
            save_untouch_nii(mu_perm_v_haufe_nii,[proj.path.haufe.fmri_ex_gm_mdl,'mu_perm_haufe_v_N=',num2str(Nperm),'_05.nii']);

        end

        % ----------------------------------------
        % Save out: mean encoding of permstrap sign. (p<0.01) group
        % gray-matter voxels
        if(numel(sig_ids_01_v)>0)

            % mu_perm_v_haufe_nii = build_nii_from_gm_mask(grp_haufe_v(sig_ids_01_v,1),gm_nii,sig_ids_01_v);
            % save_nii(mu_perm_v_haufe_nii,[proj.path.haufe.fmri_ex_gm_mdl,'mu_perm_haufe_v_N=',num2str(Nperm),'_01.nii']);

            mu_perm_v_haufe_nii = build_nii_from_gm_mask(grp_haufe_v_tstat(sig_ids_01_v,1),gm_nii,sig_ids_01_v);
            save_untouch_nii(mu_perm_v_haufe_nii,[proj.path.haufe.fmri_ex_gm_mdl,'mu_perm_haufe_v_N=',num2str(Nperm),'_01.nii']);

        end

        % ----------------------------------------
        % Save out: mean encoding of permstrap sign. (p<0.001) group gray-matter voxels
        if(numel(sig_ids_001_v)>0)
            
            % mu_perm_v_haufe_nii = build_nii_from_gm_mask(grp_haufe_v(sig_ids_001_v,1),gm_nii,sig_ids_001_v);
            %  save_nii(mu_perm_v_haufe_nii,[proj.path.haufe.fmri_ex_gm_mdl,'mu_perm_haufe_v_N=',num2str(Nperm),'_001.nii']);

            mu_perm_v_haufe_nii = build_nii_from_gm_mask(grp_haufe_v_tstat(sig_ids_001_v,1),gm_nii,sig_ids_001_v);
            save_untouch_nii(mu_perm_v_haufe_nii,[proj.path.haufe.fmri_ex_gm_mdl,'mu_perm_haufe_v_N=',num2str(Nperm),'_001.nii']);

        end

        %% ----------------------------------------
        %% Arousal
        sig_ids_05_a = [];
        sig_ids_01_a = [];
        sig_ids_001_a = [];
        
        for k=1:numel(row_ids_a)
            
            % Count extrem random samples
            Next = 0;
            if(grp_haufe_a(row_ids_a(k),1)>0)
                Next = numel(find(grp_haufe_a(row_ids_a(k),2:i)>grp_haufe_a(row_ids_a(k),1)));
            else
                Next = numel(find(grp_haufe_a(row_ids_a(k),2:i)<grp_haufe_a(row_ids_a(k),1)));
            end
            
            % Do 2-sided tests
            if(Next<round((alpha05/2)*i))
                sig_ids_05_a = [sig_ids_05_a,row_ids_a(k)];
            end
            
            if(Next<round((alpha01/2)*i))
                sig_ids_01_a = [sig_ids_01_a,row_ids_a(k)];
            end
            
            if(Next<round((alpha001/2)*i))
                sig_ids_001_a = [sig_ids_001_a,row_ids_a(k)];
            end
            
        end

        % ----------------------------------------
        % Save out: mean encoding of group gray-matter voxels
        if(numel(row_ids_a)>0)

            % mu_a_haufe_nii = build_nii_from_gm_mask(grp_haufe_a(row_ids_a,1),gm_nii,row_ids_a);
            % save_untouch_nii(mu_a_haufe_nii,[proj.path.haufe.fmri_ex_gm_mdl,'mu_haufe_a_N=',num2str(Nperm),'.nii']);

            mu_a_haufe_nii = build_nii_from_gm_mask(grp_haufe_a_tstat(row_ids_a,1),gm_nii,row_ids_a);
            save_untouch_nii(mu_a_haufe_nii,[proj.path.haufe.fmri_ex_gm_mdl,'mu_haufe_a_N=',num2str(Nperm),'.nii']);


        end

        % ----------------------------------------
        % Save out: mean encoding of permstrap sign. (p<0.05) group
        % gray-matter voxels
        if(numel(sig_ids_05_a)>0)

            % mu_perm_a_haufe_nii = build_nii_from_gm_mask(grp_haufe_a(sig_ids_05_a,1),gm_nii,sig_ids_05_a);
            % save_untouch_nii(mu_perm_a_haufe_nii,[proj.path.haufe.fmri_ex_gm_mdl,'mu_perm_haufe_a_N=',num2str(Nperm),'_05.nii']);

            mu_perm_a_haufe_nii = build_nii_from_gm_mask(grp_haufe_a_tstat(sig_ids_05_a,1),gm_nii,sig_ids_05_a);
            save_untouch_nii(mu_perm_a_haufe_nii,[proj.path.haufe.fmri_ex_gm_mdl,'mu_perm_haufe_a_N=',num2str(Nperm),'_05.nii']);

        end

        % ----------------------------------------
        % Save out: mean encoding of permstrap sign. (p<0.01) group
        % gray-matter voxels
        if(numel(sig_ids_01_a)>0)

            % mu_perm_a_haufe_nii = build_nii_from_gm_mask(grp_haufe_a(sig_ids_01_a,1),gm_nii,sig_ids_01_a);
            % save_untouch_nii(mu_perm_a_haufe_nii,[proj.path.haufe.fmri_ex_gm_mdl,'mu_perm_haufe_a_N=',num2str(Nperm),'_01.nii']);

            mu_perm_a_haufe_nii = build_nii_from_gm_mask(grp_haufe_a_tstat(sig_ids_01_a,1),gm_nii,sig_ids_01_a);
            save_untouch_nii(mu_perm_a_haufe_nii,[proj.path.haufe.fmri_ex_gm_mdl,'mu_perm_haufe_a_N=',num2str(Nperm),'_01.nii']);

        end

        % ----------------------------------------
        % Save out: mean encoding of permstrap sign. (p<0.001) group gray-matter voxels
        if(numel(sig_ids_001_a)>0)

            % mu_perm_a_haufe_nii = build_nii_from_gm_mask(grp_haufe_a(sig_ids_001_a,1),gm_nii,sig_ids_001_a);
            % save_untouch_nii(mu_perm_a_haufe_nii,[proj.path.haufe.fmri_ex_gm_mdl,'mu_perm_haufe_a_N=',num2str(Nperm),'_001.nii']);

            mu_perm_a_haufe_nii = build_nii_from_gm_mask(grp_haufe_a_tstat(sig_ids_001_a,1),gm_nii,sig_ids_001_a);
            save_untouch_nii(mu_perm_a_haufe_nii,[proj.path.haufe.fmri_ex_gm_mdl,'mu_perm_haufe_a_N=',num2str(Nperm),'_001.nii']);

        end

    end

end
