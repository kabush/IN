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
logger(['Comparing IN Cog Control Model Predictions       '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.ctrl.in_3dlme_v]);
    eval(['! rm -rf ',proj.path.ctrl.in_3dlme_v]);
    disp(['Creating ',proj.path.ctrl.in_3dlme_v]);
    eval(['! mkdir ',proj.path.ctrl.in_3dlme_v]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% Storage for analysis
all_subjs = [];
all_path = [];
all_err = [];
all_cnf = [];
all_pel = [];
all_pro = [];
all_evc = [];
all_traj = [];

%% ----------------------------------------
%% Transform beta-series into affect series {v,a}
subj_cnt = 0;
data_cnt = 0;

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

        subj_cnt = subj_cnt+1;

        %% ----------------------------------------
        %% Load computational models

        %error
        load([proj.path.ctrl.in_err_mdl,subj_study,'_',name,'_mdls.mat']);
        err_mdls = mdls;

        %conflict
        load([proj.path.ctrl.in_cnf_mdl,subj_study,'_',name,'_mdls.mat']);
        cnf_mdls = mdls;

        %prediction err likelihood
        load([proj.path.ctrl.in_pel_mdl,subj_study,'_',name,'_mdls.mat']);
        pel_mdls = mdls;

        %prediction response outcome
        load([proj.path.ctrl.in_pro_mdl,subj_study,'_',name,'_mdls.mat']);
        pro_mdls = mdls;

        %expected value of control
        load([proj.path.ctrl.in_evc_mdl,subj_study,'_',name,'_mdls.mat']);
        evc_mdls = mdls;

        %% ----------------------------------------
        %% Load fixed effects

        % true predictors
        err = reshape(sqrt((err_mdls.v_dcmp').^2),1,prod(size(err_mdls.v_dcmp)))';
        cnf = reshape(cnf_mdls.v_dcmp',1,prod(size(cnf_mdls.v_dcmp)))';
        pel = reshape(pel_mdls.v_dcmp',1,prod(size(pel_mdls.v_dcmp)))';
        pro = reshape(pro_mdls.v_dcmp',1,prod(size(pro_mdls.v_dcmp)))';
        evc = reshape(evc_mdls.v_dcmp',1,prod(size(evc_mdls.v_dcmp)))';

        %indices (using err as template)
        indx = reshape(err_mdls.v_indx',1,prod(size(err_mdls.v_indx)));

        %order of observations
        traj_box = repmat([1:4],30,1);
        traj = reshape(traj_box',1,prod(size(traj_box)))';

        % true subject id
        subj = repmat([subj_study,'_',name],numel(indx),1);

        % gather group level information
        all_err = [all_err;zscore(err)];
        all_cnf = [all_cnf;zscore(cnf)];
        all_pel = [all_pel;zscore(pel)];
        all_pro = [all_pro;zscore(pro)];
        all_evc = [all_evc;zscore(evc)];
        all_traj = [all_traj;traj];

        % gather the path information
        for j = 1:numel(indx)
            data_cnt = data_cnt + 1;
            all_path{data_cnt} = [proj.path.betas.fmri_in_beta,subj_study,'_',name,'_lss.nii''[',num2str(indx(j)),']'''];
            all_subjs{data_cnt} = [subj_study,'_',name];

        end
    end
  
end

% Grand mean center all (already zscored within subject)
all_err = all_err-mean(all_err);
all_cnf = all_cnf-mean(all_cnf);
all_evc = all_evc-mean(all_evc);
all_pel = all_pel-mean(all_pel);
all_pro = all_pro-mean(all_pro);
all_traj = all_traj-mean(all_traj);

var_names = {'err','cnf','evc','pel','pro'};

for i=1:numel(var_names) 

    %Current performance measure
    var_name = var_names{i}; 
    disp(var_name);
    
    % Build script
    fid = fopen(['./lme_',var_name,'_script'],'w');
    fprintf(fid,'#! /bin/csh\n');
    fprintf(fid,'\n');
    fprintf(fid,['3dLME -prefix lme_',var_name,' -jobs 16   \\ ']);
    fprintf(fid,['      -resid lme_',var_name,'_resid       \\ ']);
    fprintf(fid,['      -model ''traj+',var_name,'''        \\ ']);
    fprintf(fid,['      -qVars ''traj,',var_name,'''        \\ ']);
    fprintf(fid,'       -qVarCenters ''0,0''       \\ ');
    fprintf(fid,['      -ranEff ''~1+',var_name,'''         \\ ']);
    fprintf(fid,'       -mask %s \\ ',[proj.path.mri.gm_mask,'group_gm_mask.nii']);
    fprintf(fid,'       -num_glt 3                      \\ ');
    fprintf(fid,['      -gltLabel 1 ',var_name,' -gltCode 1 ''',var_name,' :'' \\']);
    fprintf(fid,['      -gltLabel 2 traj -gltCode 2 ''traj :'' \\']);
    fprintf(fid,['      -gltLabel 3 intr -gltCode 3 ''',var_name,' : 0'' \\']);
    fprintf(fid,'       -dataTable                       \\ ');
    fprintf(fid,[' Subj traj ',var_name,'  InputFile   \\ ']);

    % Write out datatable
    Nrows = size(evalin('base',['all_',var_name]),1);
    for i = 1:(Nrows-1)
        fprintf(fid,' %s %1.3f %1.3f %s   \\',all_subjs{i},all_traj(i),evalin('base',['all_',var_name,'(i)']),all_path{i});
    end
    i=Nrows;
    fprintf(fid,' %s %1.3f %1.3f %s  \n \\',all_subjs{i},all_traj(i),evalin('base',['all_',var_name,'(i)']),all_path{i});
    fclose(fid);
    
    eval(['! chmod u+x lme_',var_name,'_script']);
    eval(['! ./lme_',var_name,'_script']);
    
    eval(['! mv lme_',var_name,'+tlrc.* ',proj.path.ctrl.in_3dlme_v]);
    eval(['! mv lme_',var_name,'_resid+tlrc.* ',proj.path.ctrl.in_3dlme_v]);
    eval(['! mv lme_',var_name,'_script ',proj.path.ctrl.in_3dlme_v]);

end

%% ----------------------------------------
%% ----------------------------------------
%% CLUSTER ANALYSIS BELOW!!!
%% ----------------------------------------
%% ----------------------------------------

% eval(['! 3dcalc -a ',proj.path.ctrl.in_3dlme_v,'lme_cmb+tlrc -b ' ...
%                     ,proj.path.ctrl.in_ica, ...
%       'clst_sng_orient_thresh_zstatd70_17_3x3x3.nii ' ...
%       '-expr ''a*b'' -prefix ',proj.path.ctrl.in_3dlme_v,'lme_cmb_dACC']);



% eval(['! 3dFWHMx -acf -input ',proj.path.ctrl.in_3dlme_v,'lme_cmb_resid+tlrc -mask ',proj.path.ctrl.in_ica,...
%      'clst_sng_orient_thresh_zstatd70_17_3x3x3.nii > '])
% 
% 
% eval(['! rm 3dFWHMx.1D']);
% eval(['! rm 3dFWHMx.1D.png']);
% 
% eval(['! 3dFWHMx -acf -input ',proj.path.ctrl.in_3dlme_v, ...
%       'lme_cmb_resid+tlrc > ',proj.path.ctrl.in_3dlme_v,'smooth_params_v.txt'])

% load([proj.path.ctrl.in_3dlme_v,'smooth_params_v.txt']);
% params = smooth_params_v(2,1:3);
% 
% eval(['! 3dClustSim -acf ',...
%       num2str(params(1)),' ',...
%       num2str(params(2)),' ',...
%       num2str(params(3)),' ',...
%       ' -nxyz 54 64 50 -athr 0.05 -mask ',...
%       proj.path.mri.gm_mask,'group_gm_mask.nii']);

