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
    disp(['Removing ',proj.path.ctrl.in_3dlme_a]);
    eval(['! rm -rf ',proj.path.ctrl.in_3dlme_a]);
    disp(['Creating ',proj.path.ctrl.in_3dlme_a]);
    eval(['! mkdir ',proj.path.ctrl.in_3dlme_a]);
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
        err = reshape(sqrt((err_mdls.a_dcmp').^2),1,prod(size(err_mdls.a_dcmp)))';
        cnf = reshape(cnf_mdls.a_dcmp',1,prod(size(cnf_mdls.a_dcmp)))';
        pel = reshape(pel_mdls.a_dcmp',1,prod(size(pel_mdls.a_dcmp)))';
        pro = reshape(pro_mdls.a_dcmp',1,prod(size(pro_mdls.a_dcmp)))';
        evc = reshape(evc_mdls.a_dcmp',1,prod(size(evc_mdls.a_dcmp)))';

        %indices (using err as template)
        indx = reshape(err_mdls.a_indx',1,prod(size(err_mdls.a_indx)));

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
all_pel = all_pel-mean(all_pel);
all_pro = all_pro-mean(all_pro);
all_evc = all_evc-mean(all_evc);

% ----------------------------------------
% write out 3dLME command scripts

% error
fid = fopen('./lme_err_script','w');
fprintf(fid,'#! /bin/csh\n');
fprintf(fid,'\n');
fprintf(fid,'3dLME -prefix lme_err -jobs 1                  \\ ');
fprintf(fid,'      -resid lme_err_resid                     \\ ');
fprintf(fid,'      -model ''fixed''                         \\ ');
fprintf(fid,'      -qVars fixed                             \\ ');
fprintf(fid,'      -qVarCenters 0                           \\ ');
fprintf(fid,'      -ranEff ''~1+fixed''                     \\ ');
fprintf(fid,'      -dataTable                               \\ ');
fprintf(fid,'      Subj       fixed    InputFile   \\ ');
for i = 1:size(all_err,1)-1
    fprintf(fid,'      %s    %1.3f      %s   \\',all_subjs{i},all_err(i),all_path{i});
end
i=size(all_err,1);
fprintf(fid,'      %s    %1.3f      %s \n',all_subjs{i},all_err(i),all_path{i});
fclose(fid);
eval(['! chmod u+x lme_err_script']);



% conflict
fid = fopen('./lme_cnf_script','w');
fprintf(fid,'#! /bin/csh\n');
fprintf(fid,'\n');
fprintf(fid,'3dLME -prefix lme_cnf -jobs 1                  \\ ');
fprintf(fid,'      -resid lme_cnf_resid                     \\ ');
fprintf(fid,'      -model ''fixed''                         \\ ');
fprintf(fid,'      -qVars fixed                             \\ ');
fprintf(fid,'      -qVarCenters 0                           \\ ');
fprintf(fid,'      -ranEff ''~1+fixed''                     \\ ');
fprintf(fid,'      -dataTable                               \\ ');
fprintf(fid,'      Subj       fixed    InputFile   \\ ');
for i = 1:size(all_cnf,1)-1
    fprintf(fid,'      %s    %1.3f      %s   \\',all_subjs{i},all_cnf(i),all_path{i});
end
i=size(all_cnf,1);
fprintf(fid,'      %s    %1.3f      %s \n',all_subjs{i},all_cnf(i),all_path{i});
fclose(fid);
eval(['! chmod u+x lme_cnf_script']);

% prediction error likelihood
fid = fopen('./lme_pel_script','w');
fprintf(fid,'#! /bin/csh\n');
fprintf(fid,'\n');
fprintf(fid,'3dLME -prefix lme_pel -jobs 1                  \\ ');
fprintf(fid,'      -resid lme_pel_resid                     \\ ');
fprintf(fid,'      -model ''fixed''                         \\ ');
fprintf(fid,'      -qVars fixed                             \\ ');
fprintf(fid,'      -qVarCenters 0                           \\ ');
fprintf(fid,'      -ranEff ''~1+fixed''                     \\ ');
fprintf(fid,'      -dataTable                               \\ ');
fprintf(fid,'      Subj       fixed    InputFile   \\ ');
for i = 1:size(all_pel,1)-1
    fprintf(fid,'      %s    %1.3f      %s   \\',all_subjs{i},all_pel(i),all_path{i});
end
i=size(all_pel,1);
fprintf(fid,'      %s    %1.3f      %s \n',all_subjs{i},all_pel(i),all_path{i});
fclose(fid);
eval(['! chmod u+x lme_pel_script']);

% prediction response outcome
fid = fopen('./lme_pro_script','w');
fprintf(fid,'#! /bin/csh\n');
fprintf(fid,'\n');
fprintf(fid,'3dLME -prefix lme_pro -jobs 1                  \\ ');
fprintf(fid,'      -resid lme_pro_resid                     \\ ');
fprintf(fid,'      -model ''fixed''                         \\ ');
fprintf(fid,'      -qVars fixed                             \\ ');
fprintf(fid,'      -qVarCenters 0                           \\ ');
fprintf(fid,'      -ranEff ''~1+fixed''                     \\ ');
fprintf(fid,'      -dataTable                               \\ ');
fprintf(fid,'      Subj       fixed    InputFile   \\ ');
for i = 1:size(all_pro,1)-1
    fprintf(fid,'      %s    %1.3f      %s   \\',all_subjs{i},all_pro(i),all_path{i});
end
i=size(all_pro,1);
fprintf(fid,'      %s    %1.3f      %s \n',all_subjs{i},all_pro(i),all_path{i});
fclose(fid);
eval(['! chmod u+x lme_pro_script']);

% expected value of control
fid = fopen('./lme_evc_script','w');
fprintf(fid,'#! /bin/csh\n');
fprintf(fid,'\n');
fprintf(fid,'3dLME -prefix lme_evc -jobs 1                  \\ ');
fprintf(fid,'      -resid lme_evc_resid                     \\ ');
fprintf(fid,'      -model ''fixed''                         \\ ');
fprintf(fid,'      -qVars fixed                             \\ ');
fprintf(fid,'      -qVarCenters 0                           \\ ');
fprintf(fid,'      -ranEff ''~1+fixed''                     \\ ');
fprintf(fid,'      -dataTable                               \\ ');
fprintf(fid,'      Subj       fixed    InputFile   \\ ');
for i = 1:size(all_evc,1)-1
    fprintf(fid,'      %s    %1.3f      %s   \\',all_subjs{i},all_evc(i),all_path{i});
end
i=size(all_evc,1);
fprintf(fid,'      %s    %1.3f      %s \n',all_subjs{i},all_evc(i),all_path{i});
fclose(fid);
eval(['! chmod u+x lme_evc_script']);

% ----------------------------------------
% run commands and store data

% error
eval(['! ./lme_err_script']);
eval(['! mv lme_err+tlrc.* ',proj.path.ctrl.in_3dlme_a]);
eval(['! mv lme_err_resid+tlrc.* ',proj.path.ctrl.in_3dlme_a]);

% conflict
eval(['! ./lme_cnf_script']);
eval(['! mv lme_cnf+tlrc.* ',proj.path.ctrl.in_3dlme_a]);
eval(['! mv lme_cnf_resid+tlrc.* ',proj.path.ctrl.in_3dlme_a]);

% prediction error likelihood
eval(['! ./lme_pel_script']);
eval(['! mv lme_pel+tlrc.* ',proj.path.ctrl.in_3dlme_a]);
eval(['! mv lme_pel_resid+tlrc.* ',proj.path.ctrl.in_3dlme_a]);

% prediction response outcome
eval(['! ./lme_pro_script']);
eval(['! mv lme_pro+tlrc.* ',proj.path.ctrl.in_3dlme_a]);
eval(['! mv lme_pro_resid+tlrc.* ',proj.path.ctrl.in_3dlme_a]);

% expected value of control
eval(['! ./lme_evc_script']);
eval(['! mv lme_evc+tlrc.* ',proj.path.ctrl.in_3dlme_a]);
eval(['! mv lme_evc_resid+tlrc.* ',proj.path.ctrl.in_3dlme_a]);

% clean up
eval(['! mv lme_*_script ',proj.path.ctrl.in_3dlme_a]);

