%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================


function [] = calc_clust_thresh(proj,affect_name)

var_names = proj.param.ctrl.ccm_names; 

%% clean-up
eval(['! rm 3dFWHMx.1D']);
eval(['! rm 3dFWHMx.1D.png']);

for i=1:numel(var_names)

    var_name = var_names{i};
    disp(var_name);

    eval(['! 3dFWHMx -acf -input ',proj.path.analysis.in_3dlme, ...
          'lme_',affect_name,'_',var_name,'_resid+tlrc > ',proj.path.analysis.in_clust_thresh, ...
          'smooth_params_',affect_name,'_',var_name,'.txt'])
    
    load([proj.path.analysis.in_clust_thresh,'smooth_params_',affect_name,'_',var_name,'.txt']);
    params = eval(['smooth_params_',affect_name,'_',var_name,'(2,1:3)']);
    
    eval(['! 3dClustSim -acf ',...
          num2str(params(1)),' ',...
          num2str(params(2)),' ',...
          num2str(params(3)),' ',...
          ' -nxyz 54 64 50 -pthr 0.001 -athr 0.05 -mask ',...
          proj.path.mri.gm_mask,'group_gm_mask.nii -prefix ',...
          proj.path.analysis.in_clust_thresh,'clust_size_',affect_name,'_',var_name]);

    %% clean-up
    eval(['! rm 3dFWHMx.1D']);
    eval(['! rm 3dFWHMx.1D.png']);

end
