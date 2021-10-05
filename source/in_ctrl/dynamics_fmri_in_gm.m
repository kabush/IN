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
    disp(['Removing ',proj.path.ctrl.in_dyn]);
    eval(['! rm -rf ',proj.path.ctrl.in_dyn]);
    disp(['Creating ',proj.path.ctrl.in_dyn]);
    eval(['! mkdir ',proj.path.ctrl.in_dyn]);
end

%% ----------------------------------------
%% Load labels;
label_id = load([proj.path.trg.in,'stim_ids.txt']); %note change

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

    data_exist = 0;
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

        %% Data is present
        data_exist = 1;

        %% Vectorize the base image
        base_img = vec_img_2d_nii(base_nii);
        base_img = reshape(base_img,brain_size(1)*brain_size(2)*brain_size(3),brain_size(4));
        
        %% Concatenate the MASKED base image
        subj_img = base_img(in_brain,:)';
        
        %% Perform quality
        qlty = check_gm_img_qlty(subj_img);
        
    catch
        logger(['   -mask or beta-series does not exist'],proj.path.logfile);
    end
    
    
    if(qlty.ok & data_exist)
        
        %% Initialize the prediction structure of this subject
        prds = struct();
        prds.v_hd = zeros(numel(label_id),1);
        prds.a_hd = zeros(numel(label_id),1);
        
        mdl_exist = 0;
        try
            %% Load SVM models
            load([proj.path.mvpa.fmri_ex_gm_mdl,subj_study,'_',name,'_v_model.mat']);
            load([proj.path.mvpa.fmri_ex_gm_mdl,subj_study,'_',name,'_a_model.mat']);
            mdl_exist = 1;
        catch
            logger(['   -model does not exist.'],proj.path.logfile);
        end
            
        if(mdl_exist==1)

            %% ----------------------------------------
            %% predict IN task using EX-based models
            for j=1:numel(label_id)
                
                %% valence
                [tst_predict,hd] = predict(v_model,subj_img(j,:));
                prds.v_hd(j) = 1./(1+exp(-hd(2)));
                
                %% arousal
                [tst_predict,hd] = predict(a_model,subj_img(j,:));
                prds.a_hd(j) = 1./(1+exp(-hd(2)));
                
            end
            
            %% ----------------------------------------
            %% decompose predicted trajectories (& derivs)
            
            %% valence
            [prds.v_dcmp,prds.v_indx] = decompose_in(proj,label_id,prds.v_hd);
            
            %% arousal
            [prds.a_dcmp,prds.a_indx] = decompose_in(proj,label_id,prds.a_hd);
            
            logger('   -success',proj.path.logfile);
            
            % debug
            figure(99)
            plot(1:7,prds.v_dcmp.h(1,:));
            hold on;
            plot(1:7,prds.v_dcmp.err(1,:));
            plot(2:6,prds.v_dcmp.dh(1,:));
            plot(3:5,prds.v_dcmp.d2h(1,:));
            hold off;              
            drawnow
            
            %% Save out prediction structure
            save([proj.path.ctrl.in_dyn,subj_study,'_',name,'_prds.mat'],'prds');

        end
            
    else
        logger('   -failed quality check',proj.path.logfile);
    end
    
end

%% clean up
close all;
