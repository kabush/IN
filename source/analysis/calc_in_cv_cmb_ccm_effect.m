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
%% find voxel ids with positive z-scores
Nccm = numel(proj.param.ctrl.ccm_names);
z_imgs = zeros(numel(mask),Nccm);
f_imgs = zeros(numel(mask),Nccm);
all_pos_ids = [];

for i=1:Nccm;
    
    ccm_name = proj.param.ctrl.ccm_names{i};
    disp(ccm_name)

    z_nii = load_untouch_nii([proj.path.analysis.in_cv_cmb_clust_thresh,'dacc_clust_map_',affect_name,'_',ccm_name,'.nii']);
    z_img = vec_img_2d_nii(z_nii);
    z_imgs(:,i) = z_img;
 
    f_nii = load_untouch_nii([proj.path.analysis.in_cv_cmb_clust_thresh,'dacc_clust_map_',affect_name,'_',ccm_name,'_fstat.nii']);
    f_img = vec_img_2d_nii(f_nii);
    f_imgs(:,i) = f_img;
   
    pos_ids = find(z_img>0);
    all_pos_ids = [all_pos_ids;pos_ids];
    numel(all_pos_ids)

end
unq_pos_ids = unique(all_pos_ids);


%% ----------------------------------------
%% find max ccm at each voxel
pos_ccms = zeros(numel(unq_pos_ids),1);
for i=1:numel(unq_pos_ids);

    % voxel id
    id = unq_pos_ids(i);

    % extract scores across ccms (w/ pos betas)
    z_scores = z_imgs(id,:);
    [vals,idcs] = sort(z_scores);
    pos_ccms(i) = idcs(Nccm);

end

%% ----------------------------------------
%% build ccm id labeling of ccms
ccm_nii = build_nii_from_gm_mask(pos_ccms,gm_nii,unq_pos_ids);
save_untouch_nii(ccm_nii,[proj.path.analysis.in_cv_cmb_ccm_effect,'ccm_',affect_name,'.nii']);


%% ----------------------------------------
%% build color labeling of ccms
col_pos_ccms = pos_ccms;
if(strcmp(affect_name,'v') ~=0)
    col_pos_ccms = 2*((pos_ccms-min(pos_ccms))/(max(pos_ccms)-min(pos_ccms)))-1+0.01;
else
    col_pos_ccms = ((pos_ccms-min(pos_ccms))/(max(pos_ccms)-min(pos_ccms)))-1+0.01;
end
ccm_nii = build_nii_from_gm_mask(col_pos_ccms,gm_nii,unq_pos_ids);
save_untouch_nii(ccm_nii,[proj.path.analysis.in_cv_cmb_ccm_effect,'col_ccm_',affect_name,'.nii']);