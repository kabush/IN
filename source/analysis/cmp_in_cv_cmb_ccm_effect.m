%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================


function [] = cmp_in_cv_cmb_ccm_effect(proj)


%% ----------------------------------------
%% Load gray matter mask 
v_ccm_nii = load_untouch_nii([proj.path.analysis.in_cv_cmb_ccm_effect,'ccm_v.nii']);
v_img = vec_img_2d_nii(v_ccm_nii);

a_ccm_nii = load_untouch_nii([proj.path.analysis.in_cv_cmb_ccm_effect,'ccm_a.nii']);
a_img = vec_img_2d_nii(a_ccm_nii);

Nccm = numel(proj.param.ctrl.ccm_names);

for i=1:Nccm;

    logger([proj.param.ctrl.ccm_names{i},': '],proj.path.logfile);

    %% Get ids of positive voxels by ccm
    v_ccm_ids = find(v_img==i);
    N_v_ccm = numel(v_ccm_ids);
    a_ccm_ids = find(a_img==i);
    N_a_ccm = numel(a_ccm_ids);

    %% Find ids of voxels that preserve ccm across affect dimension
    jnt_ccm_ids = intersect(v_ccm_ids,a_ccm_ids);
    N_jnt_ccm = numel(jnt_ccm_ids);

    frac = N_jnt_ccm/N_v_ccm; %valence is the less well-represented system

    if(isnan(frac))
        logger(['    [0 valence voxels]'],proj.path.logfile);
    else
        logger(['    ',num2str(N_v_ccm),' valence voxels'],proj.path.logfile);
        logger(['    ',num2str(N_a_ccm),' arousal voxels'],proj.path.logfile);
        logger(['    ',num2str(N_jnt_ccm),' joint voxels'],proj.path.logfile);
        % logger(['    [',num2str(frac*100),'preserved]'],proj.path.logfile);
    end

end



