%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

function [] = apply_mfc_mask(proj,affect_name)

var_names = proj.param.ctrl.ccm_names; 

for i=1:numel(var_names)

    var_name = var_names{i};
    disp(var_name);

    %% Generate boolean mask for clusterized masks
    cmd = ['! 3dcalc -a ',proj.path.analysis.in_clust_thresh, ...
          'clust_mask_',affect_name,'_',...
          var_name,'+tlrc -expr ''bool(a)'' -prefix ',...
          proj.path.code,'tmp/bool_clust_mask_',affect_name,'_',...
          var_name];

    disp(cmd)
    eval(cmd)

    %% Generate boolean mask for mfc * boolean cluster mask
    cmd = ['! 3dcalc -a ',proj.path.code,'tmp/bool_clust_mask_',...
           affect_name,'_',var_name,'+tlrc -b ',...
           proj.path.ctrl.in_ica,...
           'erode2_sng_orient_mfc_mask_3x3x3.nii.gz ' ...
           '-expr ''bool(a*b)'' -prefix ',proj.path.code,...
           'tmp/mfc_clust_mask_',affect_name,'_',var_name];

    disp(cmd)
    eval(cmd)

    %% Generate *.nii and move
    cmd = ['! 3dAFNItoNIFTI ', ...
           proj.path.code,'tmp/mfc_clust_mask_',affect_name, ...
           '_',var_name,'+tlrc'];
    disp(cmd);
    eval(cmd);

    cmd = ['! mv ./mfc_clust_mask_',affect_name,'_',var_name,...
           '.nii ',proj.path.analysis.in_clust_thresh];
    disp(cmd);
    eval(cmd);

    %% Generate mfc restricted maps
    cmd = ['! 3dcalc -a ',proj.path.code,'tmp/mfc_clust_mask_',...
           affect_name,'_',var_name,'+tlrc -b ',...
           proj.path.analysis.in_clust_thresh,...
          'clust_map_',affect_name,'_',var_name,'+tlrc ' ...
           '-expr ''a*b'' -float -prefix ',proj.path.code,...
           'tmp/mfc_clust_map_',affect_name,'_',var_name];
    disp(cmd)
    eval(cmd)

    %% move
    cmd = ['! mv ',proj.path.code,'tmp/mfc_clust_map_',...
           affect_name,'_',var_name,...
           '+tlrc.* ',proj.path.analysis.in_clust_thresh];
    disp(cmd);
    eval(cmd);

    %% clean-up
    eval(['! rm ',proj.path.code,'tmp/*']);

end
