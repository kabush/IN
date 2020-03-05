%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

function [] = calc_in_cv_cmb_fmri_3dlme(proj,affect_name)

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% Storage for analysis
all_subjs = [];
all_path = [];
all_err = [];
all_pro = [];
all_cnf = [];
all_evc = [];
all_trj = [];
all_aff = [];

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
        load([proj.path.ctrl.in_dyn,subj_study,'_',name,'_prds.mat']);

        % load errors
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
        err_indx_box = eval(['mdls.',affect_name,'_indx']); 
        err = reshape(sqrt((err_box').^2),1,prod(size(err_box)))';

        % prediction response outcome
        load([proj.path.ctrl.in_pro_opt_mdl,subj_study,'_',name,'_pro_opt_',affect_name,'.mat']);
        pro = sqrt(pro_opt.^2); %rse.

        % expected value of control (cross-validated values)
        load([proj.path.ctrl.in_evc_icv_mdl,'Q_traj_',affect_name,'.mat']); 
        load([proj.path.ctrl.in_evc_icv_mdl,'Q_cv_',affect_name,'.mat']); 
        evc_box = 0*err_box;
        for j=1:numel(subjs)
            if(Q_cv(i,j)==1)
                evc_box = evc_box + squeeze(Q_traj_cv(i,j,:,:));
            end
        end
        evc_box = evc_box/sum(Q_cv(i,:));
        evc = reshape(evc_box',1,prod(size(evc_box)))';

        % conflict (based on Q-values, 1st vs 2nd best Q-values)
        load([proj.path.ctrl.in_evc_icv_mdl,'Q_cnf_',affect_name,'.mat']); 
        evc_2nd_box = 0*err_box;
        for j=1:numel(subjs)
            if(Q_cv(i,j)==1)
                evc_2nd_box = evc_2nd_box + squeeze(Q_cnf_cv(i,j,:,:));
            end
        end
        evc_2nd_box = evc_2nd_box/sum(Q_cv(i,:));
        cnf_box = 1-abs(tanh(evc_box-evc_2nd_box));  %if values very
                                                     %close, conflict
                                                     %goes up
        cnf = reshape(cnf_box',1,prod(size(cnf_box)))';

        %% ----------------------------------------
        %% Load additional effects

        % affect
        load([proj.path.ctrl.in_dyn,subj_study,'_',name,'_prds.mat']);        
        aff_box = eval(['prds.',affect_name,'_dcmp.feel']);
        aff = reshape(aff_box',1,prod(size(aff_box)))';

        % Order of observations
        traj_box = repmat([1:4],numel(err)/4,1);
        traj = reshape(traj_box',1,prod(size(traj_box)))';

        % True subject id
        subj = repmat([subj_study,'_',name],numel(err),1);
        
        % Indices of LSS
        indx_box = err_indx_box;
        indx = reshape(indx_box',1,prod(size(indx_box)))';

        %% ----------------------------------------
        %% Scale all fixed effects

        % gather group level information (z-score within subject)
        all_aff = [all_aff;zscore(aff)];
        all_err = [all_err;zscore(err)];
        all_pro = [all_pro;zscore(pro)];
        all_evc = [all_evc;zscore(evc)];
        all_cnf = [all_cnf;zscore(cnf)];
        all_trj = [all_trj;traj];

        % gather the path information (indices select correct volume)
        for j = 1:numel(indx)

            data_cnt = data_cnt + 1;

            % Load path to beta-series
            all_path{data_cnt} = [proj.path.betas.fmri_in_beta,subj_study,'_',name,'_lss.nii''[',num2str(indx(j)),']'''];
            all_subjs{data_cnt} = [subj_study,'_',name];

            % Load subject's sex & age
            demo = readtable([proj.path.raw_data,proj.path.demo,'/',subj_study,'.csv']);
            id = find(strcmp(demo.ID,name)~=0);
            sex_code = demo.Type(id);
            age = demo.Age(id);
            all_age_cell{data_cnt} = age;
            sex='N';
            if(sex_code==1)
                all_sex_cell{data_cnt} = 1;
                sex='M';
            else
                all_sex_cell{data_cnt} = -1; 
                sex='F';
            end
            all_sex_code{data_cnt} = sex;
            
        end

    end
  
end

% Convert age to matlab array
all_age = cell2mat(all_age_cell);

%% ----------------------------------------
%% base model variables
all_aff = all_aff-mean(all_aff);
all_trj = all_trj-mean(all_trj);
all_err = all_err-mean(all_err);
all_pro = all_pro-mean(all_pro);
all_cnf = all_cnf-mean(all_cnf);
all_evc = all_evc-mean(all_evc);
all_age = all_age-mean(all_age);

% Build script
fid = fopen(['./lme_',affect_name,'_cv_cmb_script'],'w');
fprintf(fid,'#! /bin/csh\n');
fprintf(fid,'\n');
fprintf(fid,['3dLME -prefix lme_',affect_name,'_cv_cmb -jobs 16 \\']);
fprintf(fid,['      -resid lme_',affect_name,'_cv_cmb_resid \\']);

%% with just sex & age & (
fprintf(fid,['      -model ''aff*sex*age+trj*sex*age+err*sex*age+pro*sex*age+cnf*sex*age+evc*sex*age'' \\']);
fprintf(fid,['      -qVars ''aff,trj,err,pro,cnf,evc,age'' \\']);
fprintf(fid,['      -qVarCenters ''0,0,0,0,0,0,0'' \\']);

fprintf(fid,['      -ranEff ''~1'' \\ ']); 
fprintf(fid,'       -mask %s \\ ',[proj.path.mri.gm_mask,'group_gm_mask.nii']);

fprintf(fid,'       -num_glt 15 \\ ');

fprintf(fid,['      -gltLabel  1  aff  -gltCode  1  ''aff :'' \\']);
fprintf(fid,['      -gltLabel  2  trj  -gltCode  2  ''trj :'' \\']);
fprintf(fid,['      -gltLabel  3  err  -gltCode  3  ''err :'' \\']);
fprintf(fid,['      -gltLabel  4  pro  -gltCode  4  ''pro :'' \\']);
fprintf(fid,['      -gltLabel  5  cnf  -gltCode  5  ''cnf :'' \\']);
fprintf(fid,['      -gltLabel  6  evc  -gltCode  6  ''evc :'' \\']);
fprintf(fid,['      -gltLabel  7  age  -gltCode  7  ''age :'' \\']);
fprintf(fid,['      -gltLabel  8  sex  -gltCode  8  ''sex : 1*M -1*F'' \\']);
fprintf(fid,['      -gltLabel  9  yint -gltCode  9  ''aff : 0'' \\']);
fprintf(fid,['      -gltLabel 10  aff_sex  -gltCode 10  ''sex : 1*M -1*F aff :'' \\']);
fprintf(fid,['      -gltLabel 11  trj_sex  -gltCode 11  ''sex : 1*M -1*F trj :'' \\']);
fprintf(fid,['      -gltLabel 12  err_sex  -gltCode 12  ''sex : 1*M -1*F err :'' \\']);
fprintf(fid,['      -gltLabel 13  pro_sex  -gltCode 13  ''sex : 1*M -1*F pro :'' \\']);
fprintf(fid,['      -gltLabel 14  cnf_sex  -gltCode 14  ''sex : 1*M -1*F cnf :'' \\']);
fprintf(fid,['      -gltLabel 15  evc_sex  -gltCode 15  ''sex : 1*M -1*F evc :'' \\']);

fprintf(fid,'       -dataTable \\ ');

%% headers
fprintf(fid,[' Subj sex age aff trj err pro cnf evc InputFile \\ ']); ...

%% Write out datatable
Nrows = size(all_err,1); 
for i = 1:(Nrows-1)
    fprintf(fid,[' %s %s %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %s \\'],...
            all_subjs{i},...
            all_sex_code{i},...
            all_age(i),...
            all_aff(i),...
            all_trj(i),...
            all_err(i),...
            all_pro(i),...
            all_cnf(i),...
            all_evc(i),...
            all_path{i});
end

i = Nrows;
    fprintf(fid,[' %s %s %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %s \n \\'],...
            all_subjs{i},...
            all_sex_code{i},...
            all_age(i),...
            all_aff(i),...
            all_trj(i),...
            all_err(i),...
            all_pro(i),...
            all_cnf(i),...
            all_evc(i),...
            all_path{i});
fclose(fid);

% Execute the script
eval(['! chmod u+x lme_',affect_name,'_cv_cmb_script']);
eval(['! ./lme_',affect_name,'_cv_cmb_script']);

% Clean-up
eval(['! mv lme_',affect_name,'_cv_cmb+tlrc.* ',proj.path.analysis.in_cv_cmb_3dlme]);
eval(['! mv lme_',affect_name,'_cv_cmb_resid+tlrc.* ',proj.path.analysis.in_cv_cmb_3dlme]);
eval(['! mv lme_',affect_name,'_cv_cmb_script ',proj.path.analysis.in_cv_cmb_3dlme]);