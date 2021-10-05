%%========================================
%%========================================
%%
%% Keith Bush, PhD (2019)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================



function [pro_all_out,pro_all,mdl] = eval_pro_cv(proj,affect_name);
    
subjs = load_subjs(proj);

% log progress
logger('* Fitting the PRO',proj.path.logfile);

pro_all = [];
pro_all_out = []; 

for i=1:numel(subjs)
    
    % extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;
    
    % log processing of subject
    logger([subj_study,'_',name],proj.path.logfile);
    
    data_exist = 0;
    try
        
        % load dynamics
        load([proj.path.ctrl.in_dyn,subj_study,'_',name,'_prds.mat']);
        
        % data is present
        data_exist = 1;
        
    catch
        logger(['   -predictions do not exist'],proj.path.logfile);
    end
    
    if(data_exist)
        
        %% ----------------------------------------
        %% Extract state from beta-series
        
        % Load gray matter mask 
        gm_nii = load_nii([proj.path.mri.gm_mask,'group_gm_mask.nii']);
        mask = double(gm_nii.img);
        brain_size=size(mask);
        mask = reshape(mask,brain_size(1)*brain_size(2)*brain_size(3),1);
        in_brain=find(mask==1);  
        
        % Load beta-series
        base_nii = load_nii([proj.path.betas.fmri_in_beta,subj_study,'_',name,'_lss.nii']);
        brain_size = size(base_nii.img);
        
        % Vectorize the base image
        base_img = vec_img_2d_nii(base_nii);
        base_img = reshape(base_img,brain_size(1)*brain_size(2)*brain_size(3),brain_size(4))';
        
        % Load ICA masks comprising the state space (and grab
        % activations)
        ica_seq = proj.param.ctrl.ica_ids; 

        Nica = numel(ica_seq);
        states = zeros(size(base_img,1),Nica);
        
        for j=1:Nica
            
            ica_n= ica_seq(j);
            
            ica_nii = load_nii([proj.path.ctrl.in_ica, ...
                                'sng_mni_thresh_zstatd20_', ...
                                num2str(ica_n),'_3x3x3.nii.gz']);
            ica = double(ica_nii.img);
            brain_size=size(ica);
            ica = reshape(ica,brain_size(1)*brain_size(2)*brain_size(3),1);
            in_brain_this_ica = find(abs(ica)>0); % 4 is the z-score threshold
                                                  % displayed in Ray 2013
            in_brain_this_ica = intersect(in_brain,in_brain_this_ica);
            states(:,j) = mean(base_img(:,in_brain_this_ica),2);
            
        end
        
        % standard score each dimension separately
        states = zscore(states);

        %% ----------------------------------------
        %% Extract labels
        indx = eval(['prds.',affect_name,'_indx.h(:,3:(end-1))']);
        err = eval(['prds.',affect_name,'_dcmp.err(:,4:end)']);
        indx_1d = reshape(indx',1,prod(size(indx)));
        err_1d = zscore(reshape(err',1,prod(size(err)))); 

        % debug
        pro_all = [pro_all,err_1d'];
        
        %% ------------------------------------------------------------
        %% Model PRO
        mdl = fitrsvm(states(indx_1d,:),err_1d);
        [err_out] = predict(mdl,states(indx_1d,:));
        err_out_1d = err_out';
        disp(['   corr=',num2str(corr(err_1d',err_out_1d'))]);

        %% ------------------------------------------------------------
        %% Store out indices
        save([proj.path.ctrl.in_pro_opt_mdl,subj_study,'_',name,'_indx_1d.mat'],'indx_1d');

        %% ------------------------------------------------------------
        %% Store out states
        save([proj.path.ctrl.in_pro_opt_mdl,subj_study,'_',name,'_states.mat'],'states');

        %% ------------------------------------------------------------
        %% Store out errors
        save([proj.path.ctrl.in_pro_opt_mdl,subj_study,'_',name,'_err_1d.mat'],'err_1d');

        %% ------------------------------------------------------------
        %% Store out model
        save([proj.path.ctrl.in_pro_opt_mdl,subj_study,'_',name,'_pro_mdl_',affect_name,'.mat'],'mdl');
 
    end
    
end


for i=1:numel(subjs)
    
    % extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    % log processing of subject
    logger([subj_study,'_',name],proj.path.logfile);
    
    try
        
        %% ----------------------------------------
        %% Load subj states/indices
        load([proj.path.ctrl.in_pro_opt_mdl,subj_study,'_',name,'_indx_1d.mat']);
        load([proj.path.ctrl.in_pro_opt_mdl,subj_study,'_',name,'_states.mat']);
        load([proj.path.ctrl.in_pro_opt_mdl,subj_study,'_',name,'_err_1d.mat']);

        cv_ids = setdiff(1:numel(subjs),i);


        pro_cv = [];
        for j = 1:numel(cv_ids)

            cv_id = cv_ids(j);
            
            % extract CV subject info
            cv_subj_study = subjs{cv_id}.study;
            cv_name = subjs{cv_id}.name;
            logger([' *',subj_study,'_',name,', CV: ',cv_subj_study,'_',cv_name],proj.path.logfile);
            
            try
                
                %% ----------------------------------------
                %% Load CV subject model
                load([proj.path.ctrl.in_pro_opt_mdl,cv_subj_study,'_',cv_name,'_pro_mdl_',affect_name,'.mat']);

                %% predict error
                [err_out] = predict(mdl,states(indx_1d,:));
                err_out_1d = err_out';
                logger(['   CV corr=',num2str(corr(err_1d',err_out_1d'))],proj.path.logfile);

                pro_cv = [pro_cv,err_out_1d'];

            catch

                disp('***CV subject results not found');

            end

        end 

        % Store out mean CV prediction for CCM
        pro_opt = mean(pro_cv,2);
        save([proj.path.ctrl.in_pro_opt_mdl,subj_study,'_',name,'_pro_opt_',affect_name,'.mat'],'pro_opt');
        
        % validate
        pro_all_out = [pro_all_out,pro_opt];

    catch
        disp('***Subject result not found');
    end

end 

toc