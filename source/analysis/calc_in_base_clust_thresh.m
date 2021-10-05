%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================


function [] = calc_in_base_clust_thresh(proj,affect_name)

%% clean-up
eval(['! rm 3dFWHMx.1D']);
eval(['! rm 3dFWHMx.1D.png']);

eval(['! 3dFWHMx -acf -input ',proj.path.analysis.in_base_3dlme, ...
      'lme_',affect_name,'_base_resid+tlrc > ',proj.path.analysis.in_base_clust_thresh, ...
      'smooth_params_',affect_name,'_base.txt'])

load([proj.path.analysis.in_base_clust_thresh,'smooth_params_',affect_name,'_base.txt']);
params = eval(['smooth_params_',affect_name,'_base (2,1:3)']);

eval(['! 3dClustSim -acf ',...
      num2str(params(1)),' ',...
      num2str(params(2)),' ',...
      num2str(params(3)),' ',...
      ' -nxyz 54 64 50 -pthr 0.001 -athr 0.05 -mask ',...
      proj.path.mri.gm_mask,'group_gm_mask.nii -prefix ',...
      proj.path.analysis.in_base_clust_thresh,'clust_size_',affect_name,'_base']);

%% clean-up
eval(['! rm 3dFWHMx.1D']);
eval(['! rm 3dFWHMx.1D.png']);

end
