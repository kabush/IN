%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

function [] = apply_in_cv_cmb_dacc_mask_2fstat(proj,affect_name)

var_names = proj.param.ctrl.ccm_f_names; 
var_ids = proj.param.ctrl.ccm_f_ids;

for i=1:numel(var_names)

    name = var_names{i};
    id = var_ids{i};

    %% Generate boolean mask for clusterized masks
    cmd = ['! 3dcalc -a ',proj.path.analysis.in_cv_cmb_clust_thresh, ...
           'clust_mask_',affect_name,'_',name,'_fstat+tlrc -expr ''bool(a)'' -prefix ',...
           proj.path.code,'tmp/bool_clust_mask_',affect_name,'_',name,'_fstat'];
    
    disp(cmd)
    eval(cmd)
    
    %% Generate boolean mask for dacc * boolean cluster mask
    cmd = ['! 3dcalc -a ',proj.path.code,'tmp/bool_clust_mask_',...
           affect_name,'_',name,'_fstat+tlrc -b ',...
           proj.path.ctrl.in_ica,...
           'mfc_restrict_sng_mni_dacc_mask_3x3x3.nii.gz ' ...
           '-expr ''bool(a*b)'' -prefix ',proj.path.code,...
           'tmp/dacc_clust_mask_',affect_name,'_',name,'_fstat'];
    
    disp(cmd)
    eval(cmd)
    
    %% Generate *.nii and move
    cmd = ['! 3dAFNItoNIFTI ', ...
           proj.path.code,'tmp/dacc_clust_mask_',affect_name, ...
           '_',name,'_fstat+tlrc'];
    disp(cmd);
    eval(cmd);
    
    cmd = ['! mv ./dacc_clust_mask_',affect_name,'_',name,'_fstat.nii ',...
           proj.path.analysis.in_cv_cmb_clust_thresh];
    disp(cmd);
    eval(cmd);
    
    %% Generate dacc restricted maps
    cmd = ['! 3dcalc -a ',proj.path.code,'tmp/dacc_clust_mask_',...
           affect_name,'_',name,'_fstat+tlrc -b ',...
           proj.path.analysis.in_cv_cmb_clust_thresh,...
           'clust_map_',affect_name,'_',name,'_fstat+tlrc ' ...
           '-expr ''a*b'' -float -prefix ',proj.path.code,...
           'tmp/dacc_clust_map_',affect_name,'_',name,'_fstat'];
    disp(cmd)
    eval(cmd)

    %% Generate *.nii and move
    cmd = ['! 3dAFNItoNIFTI ', ...
           proj.path.code,'tmp/dacc_clust_map_',affect_name, ...
           '_',name,'_fstat+tlrc'];
    disp(cmd);
    eval(cmd);
    
    cmd = ['! mv ./dacc_clust_map_',affect_name,'_',name,'_fstat.nii ',...
           proj.path.analysis.in_cv_cmb_clust_thresh];
    disp(cmd);
    eval(cmd);

    %% clean-up
    eval(['! rm ',proj.path.code,'tmp/*']);
    
end
