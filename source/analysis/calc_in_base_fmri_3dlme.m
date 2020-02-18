%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

function [] = calc_in_base_fmri_3dlme(proj,affect_name)

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% Storage for analysis
all_subjs = [];
all_path = [];
all_aff = [];
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

        % error (using ONLY for sizing)
        load([proj.path.ctrl.in_err_mdl,subj_study,'_',name,'_mdls.mat']);
        err_box = eval(['mdls.',affect_name,'_dcmp']);
        err_indx_box = eval(['mdls.',affect_name,'_indx']); 
        err = reshape(sqrt((err_box').^2),1,prod(size(err_box)))';

        %% ----------------------------------------
        %% Load additional effects

        % Affect
        load([proj.path.ctrl.in_dyn,subj_study,'_',name, ...
              '_prds.mat']);
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
        all_traj = [all_traj;traj];

        % gather the path information (indices select correct volume)
        for j = 1:numel(indx)

            data_cnt = data_cnt + 1;
            all_path{data_cnt} = [proj.path.betas.fmri_in_beta,subj_study,'_',name,'_lss.nii''[',num2str(indx(j)),']'''];
            all_subjs{data_cnt} = [subj_study,'_',name];

            % Load subject's sex & age
            demo = readtable([proj.path.raw_data,proj.path.demo,'/',subj_study,'.csv']);
            id = find(strcmp(demo.ID,name)~=0);
            sex_code = demo.Type(id);
            age = demo.Age(id);
            all_ages_cell{data_cnt} = age;
            sex='N';
            if(sex_code==1)
                sex='M';
            else
                sex='F';
            end
            all_sex{data_cnt} = sex;

        end

    end
  
end

all_age = cell2mat(all_ages_cell);

% Grand mean center all (already zscored within subject)
all_aff = all_aff-mean(all_aff);
all_traj = all_traj-mean(all_traj);
all_age = all_age-mean(all_age);

% Build script
fid = fopen(['./lme_',affect_name,'_base_script'],'w');
fprintf(fid,'#! /bin/csh\n');
fprintf(fid,'\n');
fprintf(fid,['3dLME -prefix lme_',affect_name,'_base -jobs 16 \\ ']);
fprintf(fid,['      -resid lme_',affect_name,'_base_resid \\ ']);
fprintf(fid,['      -model ''aff*sex*age+traj*sex*age'' \\']);
fprintf(fid,['      -qVars ''aff,traj,age'' \\ ']);
fprintf(fid,['      -qVarCenters ''0,0,0'' \\ ']);
fprintf(fid,['      -ranEff ''~1'' \\ ']); 
fprintf(fid,'       -mask %s \\ ',[proj.path.mri.gm_mask,'group_gm_mask.nii']);
fprintf(fid,'       -num_glt 7                      \\ ');
fprintf(fid,['      -gltLabel 1  aff   -gltCode  1  ''aff :'' \\']);
fprintf(fid,['      -gltLabel 2  traj  -gltCode  2  ''traj :'' \\']);
fprintf(fid,['      -gltLabel 3  y_int -gltCode  3  ''aff : 0'' \\']);
fprintf(fid,['      -gltLabel 4  age   -gltCode  4  ''age :'' \\']);
fprintf(fid,['      -gltLabel 5  sex   -gltCode  5  ''sex : 1*M -1*F'' \\']);
fprintf(fid,['      -gltLabel 6  aff_sex  -gltCode 6 ''sex : 1*M -1*F aff : '' \\']);
fprintf(fid,['      -gltLabel 7  traj_sex -gltCode 7 ''sex : 1*M -1*F traj : '' \\']);

fprintf(fid,'       -dataTable                       \\ ');

fprintf(fid,[' Subj sex age aff traj InputFile   \\ ']);

% Write out datatable
Nrows = size(all_traj,1); 
for i = 1:(Nrows-1)
    fprintf(fid,' %s %s %d %1.3f %1.3f %s   \\',...
            all_subjs{i},...
            all_sex{i},...
            all_age(i),...
            all_aff(i),...
            all_traj(i),...
            all_path{i});
end
i=Nrows;
    fprintf(fid,' %s %s %d %1.3f %1.3f %s \n  \\',...
            all_subjs{i},...
            all_sex{i},...
            all_age(i),...
            all_aff(i),...
            all_traj(i),...
            all_path{i});
fclose(fid);

% Execute the script
eval(['! chmod u+x lme_',affect_name,'_base_script']);
eval(['! ./lme_',affect_name,'_base_script']);

% Clean-up
eval(['! mv lme_',affect_name,'_base+tlrc.* ',proj.path.analysis.in_base_3dlme]);
eval(['! mv lme_',affect_name,'_base_resid+tlrc.* ',proj.path.analysis.in_base_3dlme]);
eval(['! mv lme_',affect_name,'_base_script ',proj.path.analysis.in_base_3dlme]);