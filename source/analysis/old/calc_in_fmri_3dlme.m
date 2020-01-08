%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

function [] = calc_in_fmri_3dlme(proj,affect_name)

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% Storage for analysis
all_subjs = [];
all_path = [];
all_err = [];
all_cnf = [];
all_cnf_alt = [];
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

        % error
        load([proj.path.ctrl.in_err_mdl,subj_study,'_',name,'_mdls.mat']);
        err_box = eval(['mdls.',affect_name,'_dcmp']);
        indx_box = eval(['mdls.',affect_name,'_indx']);

        % conflict
        load([proj.path.ctrl.in_cnf_mdl,subj_study,'_',name,'_mdls.mat']);
        cnf_box = eval(['mdls.',affect_name,'_dcmp']);

        % conflict
        load([proj.path.ctrl.in_cnf_alt_mdl,subj_study,'_',name,'_mdls.mat']);
        cnf_alt_box = eval(['mdls.',affect_name,'_dcmp']);

        % prediction err likelihood
        load([proj.path.ctrl.in_pel_mdl,subj_study,'_',name,'_mdls.mat']);
        pel_box = eval(['mdls.',affect_name,'_dcmp']);

        % prediction response outcome
        load([proj.path.ctrl.in_pro_mdl,subj_study,'_',name,'_mdls.mat']);
        pro_box = eval(['mdls.',affect_name,'_dcmp']);

        % expected value of control
        load([proj.path.ctrl.in_evc_mdl,subj_study,'_',name,'_mdls.mat']);
        evc_box = eval(['mdls.',affect_name,'_dcmp']);

        %% ----------------------------------------
        %% Load fixed effects

        % true predictors
        err = reshape(sqrt((err_box').^2),1,prod(size(err_box)))';
        cnf = reshape(cnf_box',1,prod(size(cnf_box)))';
        cnf_alt = reshape(cnf_alt_box',1,prod(size(cnf_alt_box)))';
        pel = reshape(pel_box',1,prod(size(pel_box)))';
        pro = reshape(pro_box',1,prod(size(pro_box)))';
        evc = reshape(evc_box',1,prod(size(evc_box)))';

        %indices (using err as template)
        %%%%% indx = reshape(err_mdls.v_indx',1,prod(size(err_mdls.v_indx)));
        indx = reshape(indx_box',1,prod(size(indx_box)))';

        %order of observations
        traj_box = repmat([1:4],30,1);
        traj = reshape(traj_box',1,prod(size(traj_box)))';

        % true subject id
        subj = repmat([subj_study,'_',name],numel(indx),1);

        % gather group level information
        all_err = [all_err;zscore(err)];
        all_cnf = [all_cnf;zscore(cnf)];
        all_cnf_alt = [all_cnf_alt;zscore(cnf_alt)];
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
all_cnf_alt = all_cnf_alt-mean(all_cnf_alt);
all_evc = all_evc-mean(all_evc);
all_pel = all_pel-mean(all_pel);
all_pro = all_pro-mean(all_pro);
all_traj = all_traj-mean(all_traj);

var_names = proj.param.ctrl.ccm_names;

for i=1:numel(var_names) 

    %Current performance measure
    var_name = var_names{i}; 
    disp(var_name);
    
    % Build script
    fid = fopen(['./lme_',affect_name,'_',var_name,'_script'],'w');
    fprintf(fid,'#! /bin/csh\n');
    fprintf(fid,'\n');
    fprintf(fid,['3dLME -prefix lme_',affect_name,'_',var_name,' -jobs 16   \\ ']);
    fprintf(fid,['      -resid lme_',affect_name,'_',var_name,'_resid       \\ ']);
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
    Nrows = size(eval(['all_',var_name]),1);
    for i = 1:(Nrows-1)
        fprintf(fid,' %s %1.3f %1.3f %s   \\',all_subjs{i},all_traj(i),eval(['all_',var_name,'(i)']),all_path{i});
    end
    i=Nrows;
    fprintf(fid,' %s %1.3f %1.3f %s  \n \\',all_subjs{i},all_traj(i),eval(['all_',var_name,'(i)']),all_path{i});
    fclose(fid);
    
    eval(['! chmod u+x lme_',affect_name,'_',var_name,'_script']);
    eval(['! ./lme_',affect_name,'_',var_name,'_script']);
    
    eval(['! mv lme_',affect_name,'_',var_name,'+tlrc.* ',proj.path.analysis.in_3dlme]);
    eval(['! mv lme_',affect_name,'_',var_name,'_resid+tlrc.* ',proj.path.analysis.in_3dlme]);
    eval(['! mv lme_',affect_name,'_',var_name,'_script ',proj.path.analysis.in_3dlme]);

end

