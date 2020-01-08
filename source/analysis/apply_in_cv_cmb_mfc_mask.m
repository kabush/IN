%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

function [] = apply_in_cv_cmb_mfc_mask(proj,affect_name)

%% Generate boolean mask for clusterized masks
cmd = ['! 3dcalc -a ',proj.path.analysis.in_cv_cmb_clust_thresh, ...
       'clust_mask_',affect_name,'_cv_cmb+tlrc -expr ''bool(a)'' -prefix ',...
       proj.path.code,'tmp/bool_clust_mask_',affect_name,'_cv_cmb'];

disp(cmd)
eval(cmd)

%% Generate boolean mask for mfc * boolean cluster mask
cmd = ['! 3dcalc -a ',proj.path.code,'tmp/bool_clust_mask_',...
       affect_name,'_cv_cmb+tlrc -b ',...
       proj.path.ctrl.in_ica,...
       'erode2_sng_orient_mfc_mask_3x3x3.nii.gz ' ...
       '-expr ''bool(a*b)'' -prefix ',proj.path.code,...
       'tmp/mfc_clust_mask_',affect_name,'_cv_cmb'];

disp(cmd)
eval(cmd)

%% Generate *.nii and move
cmd = ['! 3dAFNItoNIFTI ', ...
       proj.path.code,'tmp/mfc_clust_mask_',affect_name, ...
       '_cv_cmb+tlrc'];
disp(cmd);
eval(cmd);

cmd = ['! mv ./mfc_clust_mask_',affect_name,'_cv_cmb.nii ',...
       proj.path.analysis.in_cv_cmb_clust_thresh];
disp(cmd);
eval(cmd);

%% Generate mfc restricted maps
cmd = ['! 3dcalc -a ',proj.path.code,'tmp/mfc_clust_mask_',...
       affect_name,'_cv_cmb+tlrc -b ',...
       proj.path.analysis.in_cv_cmb_clust_thresh,...
       'clust_map_',affect_name,'_cv_cmb+tlrc ' ...
       '-expr ''a*b'' -float -prefix ',proj.path.code,...
       'tmp/mfc_clust_map_',affect_name,'_cv_cmb'];
disp(cmd)
eval(cmd)

%% move
cmd = ['! mv ',proj.path.code,'tmp/mfc_clust_map_',...
       affect_name,'_cv_cmb+tlrc.* ',...
       proj.path.analysis.in_cv_cmb_clust_thresh];
disp(cmd);
eval(cmd);

%% clean-up
eval(['! rm ',proj.path.code,'tmp/*']);

