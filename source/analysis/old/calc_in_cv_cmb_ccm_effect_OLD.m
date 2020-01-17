%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================


function [] = calc_in_cv_cmb_ccm_effect(proj,affect_name)

%% ----------------------------------------
%% Load gray matter mask 
gm_nii = load_untouch_nii([proj.path.mri.gm_mask,'group_gm_mask.nii']);
mask = vec_img_2d_nii(gm_nii);
in_brain=find(mask==1);

%% ----------------------------------------
%% Load F-statistic Maps
err_fstat_nii = load_untouch_nii([proj.path.analysis.in_cv_cmb_clust_thresh,'mfc_clust_map_',affect_name,'_err_fstat.nii']);
err_fstat_img = vec_img_2d_nii(err_fstat_nii);
evc_fstat_nii = load_untouch_nii([proj.path.analysis.in_cv_cmb_clust_thresh,'mfc_clust_map_',affect_name,'_evc_fstat.nii']);
evc_fstat_img = vec_img_2d_nii(evc_fstat_nii);

%% ------------------------------------------------------------
%% Queston #1: Compare model fit effects between error and EVC

logger(['Operationss to break up F-statistics by CCM'],proj.path.logfile);
err_f_ids = find(err_fstat_img(in_brain)>0);
evc_f_ids = find(evc_fstat_img(in_brain)>0);
jnt_f_ids = intersect(err_f_ids,evc_f_ids);
evc_gtr_err_f_ids = find(evc_fstat_img(in_brain(jnt_f_ids))>err_fstat_img(in_brain(jnt_f_ids)));;

%%ERR Alone
err_alone_f_ids = setdiff(err_f_ids,jnt_f_ids);
logger(['Num. sign. Error voxels: ',num2str(numel(err_alone_f_ids))],proj.path.logfile);

%%EVC Alone
evc_alone_f_ids = setdiff(evc_f_ids,jnt_f_ids);
logger(['Num. sign. EVC voxels: ',num2str(numel(evc_alone_f_ids))],proj.path.logfile);

%%EVC > ERR
evc_gtr_err_jnt_f_ids = jnt_f_ids(evc_gtr_err_f_ids);
logger(['Num. sign. voxels where f(EVC)>f(ERR): ',num2str(numel(evc_gtr_err_jnt_f_ids))],proj.path.logfile);

%%EVC < ERR
evc_lss_err_jnt_f_ids = setdiff(jnt_f_ids,evc_gtr_err_jnt_f_ids);
logger(['Num. sign. voxels where f(EVC)<f(ERR): ',num2str(numel(evc_lss_err_jnt_f_ids))],proj.path.logfile);

%% Assign winners to voxels for visualization
prt_fstat_img = 0*evc_fstat_img;
prt_fstat_img(in_brain(err_alone_f_ids)) = 1 ;
prt_fstat_img(in_brain(evc_alone_f_ids)) =  2;
prt_fstat_img(in_brain(evc_gtr_err_jnt_f_ids)) = 2;
prt_fstat_img(in_brain(evc_lss_err_jnt_f_ids)) = 1;%4
prt_fstat_nii = build_nii_from_gm_mask(prt_fstat_img(in_brain),gm_nii,in_brain);
save_untouch_nii(prt_fstat_nii,[proj.path.analysis.in_cv_cmb_ccm_effect,'prt_fstat_',affect_name,'.nii']);

%% ------------------------------------------------------------
%% Queston #2: Compare model fit effects (F-stat) between error and EVC
%% when restricted to positive and negative effects (Z-stat)

%% ----------------------------------------
%% Load Z-statistic Maps
evc_zstat_nii = load_untouch_nii([proj.path.analysis.in_cv_cmb_clust_thresh,'mfc_clust_map_',affect_name,'_evc.nii']);
evc_zstat_img = vec_img_2d_nii(evc_zstat_nii);
err_zstat_nii = load_untouch_nii([proj.path.analysis.in_cv_cmb_clust_thresh,'mfc_clust_map_',affect_name,'_err.nii']);
err_zstat_img = vec_img_2d_nii(err_zstat_nii);

%% Operationss to break up Z-statistics by CCM
logger(['Operationss to break up Z-statistics by CCM'],proj.path.logfile);

%% ERR vs EVC: Positive effects
logger(['***Positive effects***'],proj.path.logfile);
err_pos_z_ids = find(err_zstat_img(in_brain)>0);
logger(['Num. sign. (+) Error voxels: ',num2str(numel(err_pos_z_ids))],proj.path.logfile);
evc_pos_z_ids = find(evc_zstat_img(in_brain)>0);
logger(['Num. sign. (+) EVC voxels: ',num2str(numel(evc_pos_z_ids))],proj.path.logfile);
jnt_pos_z_ids = intersect(err_pos_z_ids,evc_pos_z_ids);
logger(['Num. sign. (+) Joint voxels: ',num2str(numel(jnt_pos_z_ids))],proj.path.logfile);
evc_gtr_jnt_pos_f_ids = find(evc_fstat_img(in_brain(jnt_pos_z_ids))>err_fstat_img(in_brain(jnt_pos_z_ids)));
logger(['Num. sign. (+) Joint voxels where f(EVC) > f(ERR): ',num2str(numel(evc_gtr_jnt_pos_f_ids))],proj.path.logfile);
evc_lss_jnt_pos_f_ids = find(evc_fstat_img(in_brain(jnt_pos_z_ids))<err_fstat_img(in_brain(jnt_pos_z_ids)));
logger(['Num. sign. (+) Joint voxels where f(EVC) < f(ERR): ',num2str(numel(evc_lss_jnt_pos_f_ids))],proj.path.logfile);

logger(['***EVC better than chance?'],proj.path.logfile);
Nevc_gtr_jnt_pos = numel(evc_gtr_jnt_pos_f_ids);
Nevc_lss_jnt_pos = numel(evc_lss_jnt_pos_f_ids);
pevc=Nevc_gtr_jnt_pos/(Nevc_gtr_jnt_pos+Nevc_lss_jnt_pos);
[phat,pci] = binofit(Nevc_gtr_jnt_pos,(Nevc_gtr_jnt_pos+Nevc_lss_jnt_pos));
logger(['  -Null distr. of ccm (if coin): [',num2str(pci(1)),', ',num2str(pci(2)),']']);
logger(['   -EVC prob=',num2str(pevc)],proj.path.logfile);

% Assign winners to voxels for visualization
pos_zstat_fstat_img = 0*evc_fstat_img;
pos_zstat_fstat_img(in_brain(err_pos_z_ids)) = 7;
pos_zstat_fstat_img(in_brain(evc_pos_z_ids)) =  26;
pos_zstat_fstat_img(in_brain(evc_gtr_jnt_pos_f_ids)) = 26;
pos_zstat_fstat_img(in_brain(evc_lss_jnt_pos_f_ids)) = 7;
pos_zstat_fstat_nii = build_nii_from_gm_mask(pos_zstat_fstat_img(in_brain),gm_nii,in_brain);
save_untouch_nii(pos_zstat_fstat_nii,[proj.path.analysis.in_cv_cmb_ccm_effect,'pos_zstat_fstat_',affect_name,'.nii']);

%% ERR vs EVC: Negative effects
logger(['***Negative effects***'],proj.path.logfile);
err_neg_z_ids = find(err_zstat_img(in_brain)<0);
logger(['Num. sign. (-) Error voxels: ',num2str(numel(err_neg_z_ids))],proj.path.logfile);
evc_neg_z_ids = find(evc_zstat_img(in_brain)<0);
logger(['Num. sign. (-) EVC voxels: ',num2str(numel(evc_neg_z_ids))],proj.path.logfile);
jnt_neg_z_ids = intersect(err_neg_z_ids,evc_neg_z_ids);
logger(['Num. sign. (-) Joint voxels: ',num2str(numel(jnt_neg_z_ids))],proj.path.logfile);
evc_gtr_jnt_neg_f_ids = find(evc_fstat_img(in_brain(jnt_neg_z_ids))>err_fstat_img(in_brain(jnt_neg_z_ids)));
logger(['Num. sign. (-) Joint voxels where f(EVC) > f(ERR): ',num2str(numel(evc_gtr_jnt_neg_f_ids))],proj.path.logfile);
evc_lss_jnt_neg_f_ids = find(evc_fstat_img(in_brain(jnt_neg_z_ids))<err_fstat_img(in_brain(jnt_neg_z_ids)));
logger(['Num. sign. (-) Joint voxels where f(EVC) < f(ERR): ',num2str(numel(evc_lss_jnt_neg_f_ids))],proj.path.logfile);

logger(['EVC better than chance?'],proj.path.logfile);
Nevc_gtr_jnt_neg = numel(evc_gtr_jnt_neg_f_ids);
Nevc_lss_jnt_neg = numel(evc_lss_jnt_neg_f_ids);
pevc=Nevc_gtr_jnt_neg/(Nevc_gtr_jnt_neg+Nevc_lss_jnt_neg);
[phat,pci] = binofit(Nevc_gtr_jnt_neg,(Nevc_gtr_jnt_neg+Nevc_lss_jnt_neg));
logger(['  -Null distr. of ccm (if coin): [',num2str(pci(1)),', ',num2str(pci(2)),']']);
logger(['   -EVC prob=',num2str(pevc)],proj.path.logfile);

% Assign winners to voxels for visualization
neg_zstat_fstat_img = 0*evc_fstat_img;
neg_zstat_fstat_img(in_brain(err_neg_z_ids)) = 7;
neg_zstat_fstat_img(in_brain(evc_neg_z_ids)) =  26;
neg_zstat_fstat_img(in_brain(evc_gtr_jnt_neg_f_ids)) = 26;
neg_zstat_fstat_img(in_brain(evc_lss_jnt_neg_f_ids)) = 7;
neg_zstat_fstat_nii = build_nii_from_gm_mask(neg_zstat_fstat_img(in_brain),gm_nii,in_brain);
save_untouch_nii(neg_zstat_fstat_nii,[proj.path.analysis.in_cv_cmb_ccm_effect,'neg_zstat_fstat_',affect_name,'.nii']);
