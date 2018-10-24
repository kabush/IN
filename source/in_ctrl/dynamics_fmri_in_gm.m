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
logger(['Computing IN affect dynamics'],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.ctrl.in_ctrl]);
    eval(['! rm -rf ',proj.path.ctrl.in_ctrl]);
    disp(['Creating ',proj.path.ctrl.in_ctrl]);
    eval(['! mkdir ',proj.path.ctrl.in_ctrl]);
end

%% ----------------------------------------
%% Load labels;
v_label = load([proj.path.trg.ex,'stim_v_labs.txt']);
a_label = load([proj.path.trg.ex,'stim_a_labs.txt']);
label_id = load([proj.path.trg.in,'stim_ids.txt']); %note change
v_score = load([proj.path.trg.ex,'stim_v_scores.txt']);
a_score = load([proj.path.trg.ex,'stim_a_scores.txt']);

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% ----------------------------------------
%% Transform beta-series into affect series {v,a}
for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    % log processing of subject
    logger([subj_study,'_',name],proj.path.logfile);
    
    try

        %% Load gray matter mask 
        gm_nii = load_nii([proj.path.mri.gm_mask,subj_study,'.',name,'.gm.nii']);
        mask = double(gm_nii.img);
        brain_size=size(mask);
        mask = reshape(mask,brain_size(1)*brain_size(2)*brain_size(3),1);
        in_brain=find(mask==1);  
        
        %% Load beta-series
        path = [proj.path.betas.fmri_in_beta,subj_study,'_',name,'_lss.nii'];
        base_nii = load_nii(path);
        brain_size = size(base_nii.img);
        
        %% Vectorize the base image
        base_img = vec_img_2d_nii(base_nii);
        base_img = reshape(base_img,brain_size(1)*brain_size(2)*brain_size(3),brain_size(4));
        
        %% Concatenate the MASKED base image
        subj_img = zscore(base_img(in_brain,:)'); %z-score voxels independently
        
        %% Concatenate all label/subj identifiers
        subj_id = [repmat(id,numel(label_id),1)];
        subj_i = [repmat(i,numel(label_id),1)];
        
        %% Perform quality
        qlty = check_gm_img_qlty(subj_img);
        
        if(qlty.ok)
            
            %% Initialize the prediction structure of this subject
            prds = struct();
            prds.v_hd = zeros(numel(label_id),1);
            prds.a_hd = zeros(numel(label_id),1);
            
            %% Load SVM models
            load([proj.path.mvpa.fmri_ex_gm_cls,subj_study,'_',name,'_v_model.mat']);
            load([proj.path.mvpa.fmri_ex_gm_cls,subj_study,'_',name,'_a_model.mat']);
            
            %% ----------------------------------------
            %% predict IN task using EX-based models
            for j=1:numel(label_id)
                
                %% valence
                [tst_predict,hd] = predict(v_model,subj_img(j,:));
                prds.v_hd(j) = hd(2);
                
                %% arousal
                [tst_predict,hd] = predict(a_model,subj_img(j,:));
                prds.a_hd(j) = hd(2);
                
            end
            
            %% ----------------------------------------
            %% decompose predicted trajectories (& derivs)
            
            %% valence
            prds.v_dcmp = decompose_in(proj,label_id,prds.v_hd);
            
            %% arousal
            prds.a_dcmp = decompose_in(proj,label_id,prds.a_hd);
            
            logger('   -success',proj.path.logfile);
            
            % % debug
            % figure(99)
            % plot(1:7,prds.v_dcmp.h(1,:));
            % hold on;
            % plot(1:7,prds.v_dcmp.err(1,:));
            % plot(2:6,prds.v_dcmp.dh(1,:));
            % plot(3:5,prds.v_dcmp.d2h(1,:));
            % hold off;              
            % drawnow

            %% Save out prediction structure
            save([proj.path.ctrl.in_ctrl,subj_study,'_',name,'_prds.mat'],'prds');
            
        else
            logger('   -failed quality check',proj.path.logfile);
        end

    catch
        logger('   -dynamics error: predictions not made',proj.path.logfile);
    end

end
