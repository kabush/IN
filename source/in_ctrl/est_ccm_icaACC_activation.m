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
logger(['Computing IN Cog Control EVC Models                  '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.ctrl.in_acc_activ]);
    eval(['! rm -rf ',proj.path.ctrl.in_acc_activ]);
    disp(['Creating ',proj.path.ctrl.in_acc_activ]);
    eval(['! mkdir ',proj.path.ctrl.in_acc_activ]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% Storage for analysis
all_b_evc = [];

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

        % load dynamics
        load([proj.path.ctrl.in_err_mdl,subj_study,'_',name,'_mdls.mat']);
        
        %% Data is present
        data_exist = 1;
        
    catch
        logger(['   -predictions do not exist'],proj.path.logfile);
    end

    if(data_exist)

        %% ----------------------------------------
        %% Find GM intersection with ACC ICA

        % Load GM mask
        gm_nii = load_nii([proj.path.mri.gm_mask,subj_study,'.',name,'.gm.nii']);
        mask = double(gm_nii.img);
        brain_size=size(mask);
        mask = reshape(mask,brain_size(1)*brain_size(2)*brain_size(3),1);
        in_brain=find(mask==1);  

        % Load ACCI ICA mask
        ica_nii = load_untouch_nii([proj.path.ctrl.in_ica,'sng_orient_thresh_zstatd70_17_3x3x3.nii.gz']);
        ica = double(ica_nii.img);
        brain_size=size(ica);
        ica = reshape(ica,brain_size(1)*brain_size(2)*brain_size(3),1);
        in_brain_this_ica = find(ica>0);

        % Intersection
        in_brain_this_ica = intersect(in_brain,in_brain_this_ica);
        disp(['  icaACC|GM voxels: ',num2str(numel(unique(in_brain_this_ica)))]);

        % Load beta-series
        base_nii = load_untouch_nii([proj.path.betas.fmri_in_beta,subj_study,'_',name,'_lss.nii']);
        brain_size = size(base_nii.img);
        
        % Vectorize the base image
        base_img = vec_img_2d_nii(base_nii);
        base_img = reshape(base_img,brain_size(1)*brain_size(2)*brain_size(3),brain_size(4))';

        %% ----------------------------------------
        %% Find ACC activation trajectory
        base_acc = mean(base_img(:,in_brain_this_ica),2);

        %% ----------------------------------------
        %% Load computational models

%         %error
%         load([proj.path.ctrl.in_err_mdl,subj_study,'_',name,'_mdls.mat']);
%         err_idx_1d = reshape(mdls.v_indx.err',1,prod(size(mdls.v_indx.err)));
%         err_v_1d = reshape(mdls.v_dcmp.err',1,prod(size(mdls.v_dcmp.err)));
%         figure(1)
%         scatter(base_acc(err_idx_1d),err_v_1d);
%         hold on;
%         [b stat] = robustfit(base_acc(err_idx_1d),err_v_1d);
%         plot(sort(base_acc(err_idx_1d)),sort(base_acc(err_idx_1d))*b(2)+b(1));
%         disp(['p=',num2str(stat.p(2))]);
% 
%         %conflict
%         load([proj.path.ctrl.in_cnf_mdl,subj_study,'_',name,'_mdls.mat']);
%         cnf_idx_1d = reshape(mdls.v_indx.cnf',1,prod(size(mdls.v_indx.cnf)));
%         cnf_v_1d = reshape(mdls.v_dcmp.cnf',1,prod(size(mdls.v_dcmp.cnf)));
%         figure(2)
%         scatter(base_acc(cnf_idx_1d),cnf_v_1d);
%         hold on;
%         [b stat] = robustfit(base_acc(cnf_idx_1d),cnf_v_1d);
%         plot(sort(base_acc(cnf_idx_1d)),sort(base_acc(cnf_idx_1d))*b(2)+b(1));
%         disp(['p=',num2str(stat.p(2))]);
% 
%         %pel
%         load([proj.path.ctrl.in_pel_mdl,subj_study,'_',name,'_mdls.mat']);
%         pel_idx_1d = reshape(mdls.v_indx.pel',1,prod(size(mdls.v_indx.pel)));
%         pel_v_1d = reshape(mdls.v_dcmp.pel',1,prod(size(mdls.v_dcmp.pel)));
%         figure(3)
%         scatter(base_acc(pel_idx_1d),pel_v_1d);
%         hold on;
%         [b stat] = robustfit(base_acc(pel_idx_1d),pel_v_1d);
%         plot(sort(base_acc(pel_idx_1d)),sort(base_acc(pel_idx_1d))*b(2)+b(1));
%         disp(['p=',num2str(stat.p(2))]);
% 
%         %pro
         load([proj.path.ctrl.in_pro_mdl,subj_study,'_',name,'_mdls.mat']);
         pro_idx_1d = reshape(mdls.v_indx.pro',1,prod(size(mdls.v_indx.pro)));
%         pro_v_1d = reshape(mdls.v_dcmp.pro',1,prod(size(mdls.v_dcmp.pro)));
%         figure(4)
%         scatter(base_acc(pro_idx_1d),pro_v_1d);
%         hold on;
%         [b stat] = robustfit(base_acc(pro_idx_1d),pro_v_1d);
%         plot(sort(base_acc(pro_idx_1d)),sort(base_acc(pro_idx_1d))*b(2)+b(1));
%         disp(['p=',num2str(stat.p(2))]);

        %evc
        load([proj.path.ctrl.in_evc_mdl,subj_study,'_',name,'_mdls.mat']);
        evc_idx_1d = pro_idx_1d; %reshape(mdls.v_indx.evc',1,evc(size(mdls.v_indx.evc)));
        evc_v_1d = reshape(mdls.v_dcmp.evc',1,prod(size(mdls.v_dcmp.evc)));
        figure(5)
        scatter(base_acc(evc_idx_1d),evc_v_1d);
        hold on;
        [b stat] = robustfit(base_acc(evc_idx_1d),evc_v_1d);
        all_b_evc = [all_b_evc,b(2)];

        plot(sort(base_acc(evc_idx_1d)),sort(base_acc(evc_idx_1d))*b(2)+b(1));
        ylim([-5,-1]);
        xlim([-1,1]);

        drawnow;
        pause(1);

        disp(['b=',num2str(b(2)),', p=',num2str(stat.p(2))]);

    end
  
end
