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
logger(['Computing IN Cog Control EVC Models                  '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.ctrl.in_evc_mdl]);
    eval(['! rm -rf ',proj.path.ctrl.in_evc_mdl]);
    disp(['Creating ',proj.path.ctrl.in_evc_mdl]);
    eval(['! mkdir ',proj.path.ctrl.in_evc_mdl]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% ----------------------------------------
%% Transform beta-series into affect series {v,a}
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
        
        %% Data is present
        data_exist = 1;
        
    catch
        logger(['   -predictions do not exist'],proj.path.logfile);
    end

    if(data_exist)

        %% Initialize the prediction structure of this subject
        mdls = struct();

        %% ----------------------------------------
        %% Build tuples for RL
        s_indx = prds.v_indx.h(:,3:(end-1));
        actions = prds.v_dcmp.h(:,4:end)-prds.v_dcmp.h(:,3:(end-1));
        rewards = -abs(actions)-abs(prds.v_dcmp.err(:,4:end));

        %% ----------------------------------------
        %% Extract state from beta-series
        % Load gray matter mask 
        %gm_nii = load_nii([proj.path.mri.gm_mask,subj_study,'.',name,'.gm.nii']);
        gm_nii = load_nii([proj.path.mri.gm_mask,'group_gm_mask.nii']);
        mask = double(gm_nii.img);
        brain_size=size(mask);
        mask = reshape(mask,brain_size(1)*brain_size(2)*brain_size(3),1);
        in_brain=find(mask==1);  

        % Load beta-series
        base_nii = load_nii([proj.path.betas.fmri_in_beta,subj_study,'_',name,'_lss.nii']);
        brain_size = size(base_nii.img);
        
        % Vectorize the base image
        base_img = vec_img_2d_nii(base_nii);
        base_img = reshape(base_img,brain_size(1)*brain_size(2)*brain_size(3),brain_size(4))';
        
        % Load ICA masks comprising the state space
        Nica = 5;

        state = zeros(size(base_img,1),Nica);

        for j=1:Nica

            ica_nii = load_nii([proj.path.ctrl.in_ica, ...
                                'orient_thresh_zstatd20_', ...
                                num2str(j),'.nii.gz']);
            ica = double(ica_nii.img);
            brain_size=size(ica);
            ica = reshape(ica,brain_size(1)*brain_size(2)*brain_size(3),1);

            in_brain_this_ica = find(abs(ica)>0); % 4 is the z-score threshold
                                                  % displayed in Ray 2013
            disp(['before: ',num2str(numel(unique(in_brain_this_ica)))]);

            in_brain_this_ica = intersect(in_brain,in_brain_this_ica);
            disp(['  after: ', ...
                  num2str(numel(unique(in_brain_this_ica)))]);

            state(:,j) = mean(base_img(:,in_brain_this_ica),2);

        end

        % standard score each dimension separately
        state = zscore(state);

        % figure(1);
        % plot(state);
        % drawnow;
        % pause(3);
        
    end
    
end
