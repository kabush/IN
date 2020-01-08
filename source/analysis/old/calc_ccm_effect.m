%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

function [] = calc_ccm_effect(proj,affect_name)

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

var_names = proj.param.ctrl.ccm_names; 

for i=1:numel(var_names)

    var_name = var_names{i};
    disp(var_name);

    %% Storage for analysis
    measures = [];
    prds = [];
    subjects = [];
    trajs = [];
    
    %% ----------------------------------------
    %% Transform beta-series into affect series {v,a}
    subj_cnt = 0;
    for j = 1:numel(subjs)
        
        %% extract subject info
        subj_study = subjs{j}.study;
        name = subjs{j}.name;
        id = subjs{j}.id;
        
        % log processing of subject
        logger([subj_study,'_',name],proj.path.logfile);
        
        data_exist = 0;
        try
            
            % load dynamics
            fpath = eval(['proj.path.ctrl.in_',var_name,'_mdl']);
            full_fpath = ([fpath,subj_study,'_',name,'_mdls.mat']);
            load(full_fpath);

            % Load mFC Mask
            mfc_nii = load_nii([proj.path.mri.gm_mask,subj_study,'.',name,'.gm.nii']);
            mfc_mask = double(mfc_nii.img);
            mfc_brain_size=size(mfc_mask);
            mfc_mask = reshape(mfc_mask,mfc_brain_size(1)*mfc_brain_size(2)*mfc_brain_size(3),1);
            mfc_in_brain=find(mfc_mask==1);  
            
            % Load mFC-masked cluster mask
            clust_nii = ...
                load_untouch_nii([proj.path.analysis ...
                                .in_clust_thresh,'mfc_clust_mask_',affect_name,'_',var_name,'.nii']); 
            clust = double(clust_nii.img);
            clust_brain_size=size(clust);
            clust = reshape(clust,clust_brain_size(1)*clust_brain_size(2)*clust_brain_size(3),1);
            clust_in_brain = find(clust==1);
            
            % Intersection
            in_brain = intersect(mfc_in_brain,clust_in_brain);
            disp(['  mFC | clust voxels: ',num2str(numel(unique(in_brain)))]);;

            % Load beta-series
            base_nii = load_untouch_nii([proj.path.betas.fmri_in_beta,subj_study,'_',name,'_lss.nii']);
            brain_size = size(base_nii.img);
            
            % Vectorize the base image
            base_img = vec_img_2d_nii(base_nii);
            base_img = reshape(base_img,brain_size(1)*brain_size(2)*brain_size(3),brain_size(4))';


            %% Data is present
            data_exist = 1;
            
        catch
            logger(['   -predictions do not exist'],proj.path.logfile);
        end
        
        if(data_exist)
            
            subj_cnt = subj_cnt+1;
            
            %% ----------------------------------------
            %% Find mFC activation trajectory
            base_mfc = mean(base_img(:,in_brain),2);
            
            %% ----------------------------------------
            %% Load computational model (named mdls)

            fpath = eval(['proj.path.ctrl.in_',var_name,'_mdl']);
            full_fpath = ([fpath,subj_study,'_',name,'_mdls.mat']);
            load(full_fpath);

            % lss indices
            if(strcmp(var_name,'evc')~=0)
                disp('EVC!');
                cmd = ['reshape(mdls.',affect_name,'_indx.evc'',1,prod(size(mdls.',affect_name,'_indx.evc)));'];
            else
                cmd = ['reshape(mdls.',affect_name,'_indx'',1,prod(size(mdls.',affect_name,'_indx)));'];
            end
            indx = eval(cmd);
            
            % true measure
            mfc = base_mfc(indx,1);
            
            % true predictors
            cmd = ['reshape(mdls.',affect_name,'_dcmp'',1,prod(size(mdls.',affect_name,'_dcmp)))'';'];
            prd = eval(cmd);
            
            % adust for error to square
            if(strcmp(var_name,'err')~=0)
                disp('ERROR!');
                prd = sqrt(prd.^2);
            end
            
            [Nkpt,Nrep] = eval(['size(mdls.',affect_name,'_dcmp)']);
            
            traj_box = repmat([1:Nrep],Nkpt,1);
            traj = reshape(traj_box',1,prod(size(traj_box)))';
            
            % true subject id
            subject = repmat(subj_cnt,numel(indx),1);
            
            % concatenate
            measures = [measures;zscore(mfc)];
            prds = [prds;zscore(prd)];
            subjects = [subjects;subject];
            trajs = [trajs;traj];
            
        end
        
    end
    
    %% ----------------------------------------
    %% Assemble Measures
    measures = double(measures-mean(measures));
    prds = double(prds-mean(prds));
    trajs = double(trajs-mean(trajs));
    subjects = double(subjects);

    size(measures)
    size(prds)
    size(trajs)
    size(subjects)
    
    %% Group GLMM fit
    tbl = table(measures,prds,trajs,subjects,'VariableNames',{'trg', ...
                        'prd','traj','subj'});;
    
    mdl_fe = fitlme(tbl,['trg ~ 1 + prd + traj']);
    mdl_re = fitlme(tbl,['trg ~ 1 + prd + traj + (prd|subj)']);

    fe_vs_re = compare(mdl_fe,mdl_re);
    mdl = mdl_fe;
    if(fe_vs_re.pValue<0.05)
        mdl = mdl_re;
        logger(' ---Random effects matter',proj.path.logfile);
    end
    
    %% ----------------------------------------
    %% compute effect size
    Rsqr = mdl.Rsquared.Ordinary;
    Fsqr = Rsqr/(1-Rsqr);
    logger(['Overall model fit: ',affect_name,'_',var_name],proj.path.logfile);
    logger(['  Rsqr=',num2str(Rsqr)],proj.path.logfile);
    logger(['  Fsqr=',num2str(Fsqr)],proj.path.logfile);
    disp(' ');

    mdl

end
