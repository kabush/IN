%%========================================
%%========================================
%%
%% Keith Bush, PhD (2019)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

function [Q_traj_cv,Q_rand_cv,act_err_cv,eval_sbj_ids] = ...
    eval_qfunc_param_cv(proj, ...
                        act_5part, ...
                        grp_act_mu, ...
                        grp_act_std, ...
                        grp_err_mu, ...
                        grp_err_std, ...
                        rand_subj_ids,gamma, ...
                        rwrd_act_f,affect_name);

% load subjs
subjs = load_subjs(proj);

% number of subjs to evaluat params
Nrand = numel(rand_subj_ids);

% log progress
logger('* Fitting the Q-functions',proj.path.logfile);

for i=1:Nrand
    
    sbj_id = rand_subj_ids(i);
    
    % extract subject info
    subj_study = subjs{sbj_id}.study;
    name = subjs{sbj_id}.name;
    id = subjs{sbj_id}.id;
    
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
        ica_seq = proj.param.ctrl.ica_ids; 
        Nica = numel(ica_seq);;
        all_states = zeros(size(base_img,1),Nica);
        
        for j=1:Nica
            
            ica_n= ica_seq(j);
            
            ica_nii = load_nii([proj.path.ctrl.in_ica, ...
                                'sng_mni_thresh_zstatd20_', ...
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
        actions_1d = zscore(reshape(actions',1,prod(size(actions))));
        errors_1d = zscore(reshape(errors',1,prod(size(errors)))); 
        
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
        %% format data for writing out structure
        %%   N   // number of tuples
        %%   X   // states
        %%   U   // actions
        %%   Xp  // next states
        %%   R   // reward (at Xp)
        %%   T   // Xp terminal state (1) or continuing (0)?
        N = numel(s_indx_1d);
        X = all_states(s_indx_1d,:)';
        Xp = all_states(sp_indx_1d,:)';
        R = rewards_1d;
        T = terminals_1d;
        U = dsc_actions_1d;
        
        % load into struct
        samples = varstostruct('N','X','U','Xp','R','T');
        
        %% ----------------------------------------
        %% configure fitted Q-iteration
        cfg = struct;
        cfg.run = 1;  % set the run flag;
        cfg.gamma = gamma; % discount factor
        cfg.U = unique(U);
        cfg.datadir = [proj.path.code,'tmp'];
        cfg.datafile = [proj.path.ctrl.in_evc_opt_mdl,subj_study,'_',name,'_result_',affect_name];
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
        
    end
    
end

%% ----------------------------------------
%% Analyze control for sign. existence
Q_traj_cv = zeros(Nrand,30,4);
Q_rand_cv = zeros(Nrand,30,4);
act_err_cv = zeros(Nrand,30,4);

sbj_i = 0;
eval_sbj_ids = [];

% log progress
logger('* CV of Fit Q-functions',proj.path.logfile);

for i=1:Nrand
    
    sbj_id = rand_subj_ids(i);
    
    % extract subject info
    subj_study = subjs{sbj_id}.study;
    name = subjs{sbj_id}.name;
    id = subjs{sbj_id}.id;

    % log processing of subject
    logger([subj_study,'_',name],proj.path.logfile);
    
    try
        
        %% ----------------------------------------
        %% Load subject (state and actions)
        load([proj.path.ctrl.in_evc_opt_mdl,subj_study,'_',name,'_result_',affect_name,'.mat']);
        this_U = Us;
        this_X = Xs';
        
        Q_traj_sbj = [];
        Q_rand_sbj = [];
        act_err_sbj = [];
        
        cv_ids = setdiff(rand_subj_ids,sbj_id);
        
        for vv = 1:numel(cv_ids)
            
            jj = cv_ids(vv);
            
            %% extract subject info
            cv_subj_study = subjs{jj}.study;
            cv_name = subjs{jj}.name;
            
            try
                
                %% ----------------------------------------
                %% Load CV subject (Q-functions)
                load([proj.path.ctrl.in_evc_opt_mdl,cv_subj_study,'_',cv_name,'_result_',affect_name,'.mat']);
                this_reg = reg;
                
                Q_traj = [];
                Q_rand = [];
                policy = [];
                opt_policy = [];
                
                parfor k=1:numel(this_U)
                    
                    u = find(this_U(k)==cfg.U);
                    q = rtenspred(this_reg{u},this_X(k,:));
                    policy = [policy,this_U(k)];
                    Q_traj = [Q_traj;q];
                    
                    %% Compute the random action as the average                            
                    Naction=numel(unique(U));
                    Q_values = [];
                    Q_run = 0;
                    
                    for m=1:Naction
                        pu = numel(find(this_U==cfg.U(m)))/numel(this_U);
                        q = rtenspred(this_reg{m},this_X(k,:));
                        Q_values = [Q_values,q];
                        Q_run = Q_run + pu*q;
                    end
                    opt_policy = [opt_policy,U(find(Q_values==max(Q_values)))];
                    Q_rand = [Q_rand;Q_run];
                    
                end
                
                Q_traj_sbj = [Q_traj_sbj; Q_traj'];
                Q_rand_sbj = [Q_rand_sbj; Q_rand'];
                act_err_sbj = [act_err_sbj; policy-opt_policy];
                
            catch
                logger([' CV: ',cv_subj_study,'_',cv_name,' results not found'],proj.path.logfile);
            end
            
        end %jj
        
        sbj_i = sbj_i+1;
        Q_traj_cv(sbj_i,:,:) = reshape(median(Q_traj_sbj),4,30)';
        Q_rand_cv(sbj_i,:,:) = reshape(median(Q_rand_sbj),4,30)';
        act_err_cv(sbj_i,:,:) = reshape(median(act_err_sbj),4,30)';                

        eval_sbj_ids = [eval_sbj_ids,sbj_id];
        
    catch
        disp('  subject result not found');
    end
    
end %i


toc