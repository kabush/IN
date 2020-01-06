%%========================================
%%========================================
%%
%% Keith Bush, PhD (2019)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

function [dcmp,indx] = eval_qfunc_cv(proj, ...
                                     subj, ...
                                     act_5part, ...
                                     grp_act_mu, ...
                                     grp_act_std, ...
                                     grp_err_mu, ...
                                     grp_err_std, ...
                                     gamma, ...
                                     rwrd_act_f,...
                                     affect_name)
    
q_traj = [];
dcmp = [];
indx = [];
    
% extract subject info
subj_study = subj.study;
name = subj.name;
id = subj.id;

% log processing of subject
logger([subj_study,'_',name],proj.path.logfile);

data_exist = 0;
try
    
    % load dynamics
    load([proj.path.ctrl.in_dyn,subj_study,'_',name,'_prds.mat']);
    
    % data is present
    data_exist = 1;
    
catch
    logger(['   -predictions do not exist'],proj.path.logfile);
end

if(data_exist)
    
    %% ----------------------------------------
    %% Extract state from beta-series
    
    % Load gray matter mask 
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
    
    % Load ICA masks comprising the state space (and grab
    % activations)
    ica_seq = [1:5]; % *** TICKET: EMOTION ICs ONLY ***
    Nica = numel(ica_seq);
    all_states = zeros(size(base_img,1),Nica);
    
    for j=1:Nica
        
        ica_n= ica_seq(j);
        
        ica_nii = load_nii([proj.path.ctrl.in_ica, ...
                            'sng_orient_thresh_zstatd20_', ...
                            num2str(ica_n),'_3x3x3.nii.gz']);
        ica = double(ica_nii.img);
        brain_size=size(ica);
        ica = reshape(ica,brain_size(1)*brain_size(2)*brain_size(3),1);
        in_brain_this_ica = find(abs(ica)>0); % 4 is the z-score threshold
                                              % displayed in Ray 2013
        in_brain_this_ica = intersect(in_brain,in_brain_this_ica);
        all_states(:,j) = mean(base_img(:,in_brain_this_ica),2);
        
    end
    
    % standard score each dimension separately
    all_states = zscore(all_states);
    
    %% ------------------------------------------------------------
    %% Calculate Q-function
    
    % Build tuples for RL
    s_indx = eval(['prds.',affect_name,'_indx.h(:,3:(end-1))']);
    sp_indx = eval(['prds.',affect_name,'_indx.h(:,4:end)']);
    pre_vals = eval(['prds.',affect_name,'_dcmp.h(:,3:(end-1))']);
    post_vals = eval(['prds.',affect_name,'_dcmp.h(:,4:end)']);
    actions = post_vals-pre_vals;
    errors = eval(['prds.',affect_name,'_dcmp.err(:,4:end)']);
    terminals = zeros(size(s_indx));terminals(:,end) = 1;
    
    % reshape to 1D
    s_indx_1d = reshape(s_indx',1,prod(size(s_indx)));
    sp_indx_1d = reshape(sp_indx',1,prod(size(sp_indx)));
    terminals_1d = reshape(terminals',1,prod(size(terminals)));
    actions_1d = zscore(reshape(actions',1,prod(size(actions))));%TICKET
    errors_1d = zscore(reshape(errors',1,prod(size(errors)))); %TICKET

    % assign action to discrete space (5 actions)
    dsc_actions_1d = 0*actions_1d;
    dsc_actions_1d(find(actions_1d>act_5part(3)))=1; 
    dsc_actions_1d(find(actions_1d>act_5part(4)))=2; 
    dsc_actions_1d(find(actions_1d<act_5part(2)))=-1; 
    dsc_actions_1d(find(actions_1d<act_5part(1)))=-2; 
    
    % Rescale actions and errors at the group level
    rwd_actions_1d = (dsc_actions_1d-grp_act_mu)/grp_act_std;
    rwd_errors_1d = (errors_1d-grp_err_mu)/grp_err_std;
    
    %% compute reward
    f_cost_act = rwrd_act_f; 
    f_cost_err = 1;
    rewards_1d = -f_cost_act*abs(rwd_actions_1d)-f_cost_err*sqrt(rwd_errors_1d.^2);
    

    %% ----------------------------------------
    %% ----------------------------------------
    %% ----------------------------------------
    for j=1:30 %%***TICKET Hardcode number of trials
        
        sbj_ids = (((j-1)*4)+1):(j*4);
        cv_ids = setdiff(1:120,sbj_ids);

        %% ----------------------------------------
        %% format data for writing out structure
        %%   N   // number of tuples
        %%   X   // states
        %%   U   // actions
        %%   Xp  // next states
        %%   R   // reward (at Xp)
        %%   T   // Xp terminal state (1) or continuing (0)?
        N = numel(s_indx_1d(cv_ids));
        X = all_states(s_indx_1d(cv_ids),:)';
        Xp = all_states(sp_indx_1d(cv_ids),:)';
        R = rewards_1d(cv_ids);
        T = terminals_1d(cv_ids);
        U = dsc_actions_1d(cv_ids);

        % for cv
        Xsbj = all_states(s_indx_1d(sbj_ids),:)';
        Usbj = dsc_actions_1d(sbj_ids);
        
        % load into struct
        samples = varstostruct('N','X','U','Xp','R','T');
        
        %% ----------------------------------------
        %% configure fitted Q-iteration
        cfg = struct;
        cfg.run = 1;  % set the run flag;
        cfg.gamma = gamma; % discount factor
        cfg.U = unique(U);
        cfg.datadir = [proj.path.code,'tmp'];
        cfg.datafile = [proj.path.ctrl.in_evc_cv_mdl,subj_study,'_',name,'_result_',affect_name,'_',num2str(j)];
        cfg.maxiter = 1000;
        cfg.regmethod = 'extratrees';
        cfg.singlereg = 0;
        cfg.trees_ntrees = 25;
        cfg.trees_nmin = 2;
        cfg.trees_k = size(X,2)+1;
        cfg.samples = samples;
        
        %% ----------------------------------------
        %% execute fitted Q-iteration
        fittedqi(cfg);
        
        %% ----------------------------------------
        %% load model and predict outlier trial
        try

                        
            %% ----------------------------------------
            %% Load CV subject (Q-functions)
            load([proj.path.ctrl.in_evc_cv_mdl,subj_study,'_',name,'_result_',affect_name,'_',num2str(j),'.mat']);

            this_U = Usbj';
            this_X = Xsbj';
            this_reg = reg;
                
            for k=1:numel(Usbj)
                    
                u = find(this_U(k)==cfg.U);
                q = rtenspred(this_reg{u},this_X(k,:));
                q_traj = [q_traj;q];

            end %% k
        
        catch
            disp('  subject result not found');
        end

    end

    %% ----------------------------------------
    %% ----------------------------------------
    %% ----------------------------------------

    %% ----------------------------------------
    %% construct models
    dcmp = q_traj; 
    dcmp = reshape(dcmp,size(s_indx,2),size(s_indx,1))';
    indx = s_indx;
    
end
