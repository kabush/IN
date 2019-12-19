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
logger(['Computing IN Cog Control EVC Optimal Parameters  '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.ctrl.in_evc_opt_mdl]);
    eval(['! rm -rf ',proj.path.ctrl.in_evc_opt_mdl]);
    disp(['Creating ',proj.path.ctrl.in_evc_opt_mdl]);
    eval(['! mkdir ',proj.path.ctrl.in_evc_opt_mdl]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% ----------------------------------------
%% discretize action set
Nsbj = 0;
act_v_all = [];
act_a_all = [];
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

        %% ******* ACTIONS ********

        %% valence
        actions = (prds.v_dcmp.h(:,4:end)-prds.v_dcmp.h(:,3:(end-1)));
        actions_1d = zscore(reshape(actions',1,prod(size(actions))));%TICKET        
        act_v_all = [act_v_all,actions_1d];

        %% arousal
        actions = (prds.a_dcmp.h(:,4:end)-prds.a_dcmp.h(:,3:(end-1)));
        actions_1d = zscore(reshape(actions',1,prod(size(actions))));%TICKET        
        act_a_all = [act_a_all,actions_1d];

        %% count subjects
        Nsbj = Nsbj + 1;

    end

end

%% Use group action to build 3 & 5 discrete action partitions
N_5size = round(numel(act_v_all)/5);

act_v_srt = sort(act_v_all);
act_a_srt = sort(act_a_all);

%% compute standard deviation
act_v_std = std(act_v_srt);

%% compute partition based on standard deviations
act_5part = [-2*act_v_std,-act_v_std,act_v_std,2*act_v_std];

%% ----------------------------------------
%% Balance action and error scales
dsc_act_v_all = [];
dsc_err_v_all = [];

for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    % log processing of subject
    logger([subj_study,'_',name],proj.path.logfile);

    data_exist = 0;
    try

        %% load dynamics
        load([proj.path.ctrl.in_dyn,subj_study,'_',name,'_prds.mat']);
        
        %% Data is present
        data_exist = 1;
        
    catch
        logger(['   -predictions do not exist'],proj.path.logfile);
    end

    if(data_exist)

        %% ******* ACTIONS ********

        %% valence
        actions = (prds.v_dcmp.h(:,4:end)-prds.v_dcmp.h(:,3:(end-1)));
        actions_1d = zscore(reshape(actions',1,prod(size(actions))));

        %% 5 discrete actions
        dsc_actions_1d = 0*actions_1d;
        dsc_actions_1d(find(actions_1d>act_5part(3)))=1; 
        dsc_actions_1d(find(actions_1d>act_5part(4)))=2; 
        dsc_actions_1d(find(actions_1d<act_5part(2)))=-1; 
        dsc_actions_1d(find(actions_1d<act_5part(1)))=-2; 
        dsc_act_v_all = [dsc_act_v_all,dsc_actions_1d];

        %% ******* ERRORS ********
        errors = (prds.v_dcmp.err(:,4:end));
        errors_1d = zscore(reshape(errors',1,prod(size(errors))));
        dsc_err_v_all = [dsc_err_v_all,errors_1d];
        
    end

end

% Compute group params
grp_act_v_mean = mean(dsc_act_v_all);
grp_act_v_std = std(dsc_act_v_all);
grp_err_v_mean = mean(dsc_err_v_all);
grp_err_v_std = std(dsc_err_v_all);

%% ----------------------------------------
%% ----------------------------------------
%% Compute Q-functions for all subjects
%% ----------------------------------------
%% ----------------------------------------

%% Meta RL Parameter search
action_5dscr_set = [1,0]; % otherwise 3 discrete actions
discount_set = [0:.1:1];
reward_act_set = [0:.2:1]; % balance between reward/action

%% Calculate set-sizes
Nact = numel(action_5dscr_set);
Ndsct = numel(discount_set);
Nfrac = numel(reward_act_set);

Q_traj_v_all = zeros(Nact,Ndsct,Nfrac,Nsbj,30,4);
Q_rand_v_all = zeros(Nact,Ndsct,Nfrac,Nsbj,30,4);
act_err_v_all = zeros(Nact,Ndsct,Nfrac,Nsbj,30,4);
Q_best_v_all = zeros(Nact,Ndsct,Nfrac,Nsbj,30,4);
Q_worst_v_all = zeros(Nact,Ndsct,Nfrac,Nsbj,30,4);

tic

%% META-LOOPS HERE
a = 1; %%only using 5 discrete actions
 
for b=6:Ndsct %% *** debug ***
    
    for c=1:Nfrac
        
        act_dscr = action_5dscr_set(a);
        gamma = discount_set(b);
        rwrd_act_f = reward_act_set(c);
        
        %% ----------------------------------------
        %% ----------------------------------------
        %% FIT Q-function to the IN task
        %% ----------------------------------------
        %% ----------------------------------------
        for i=1:numel(subjs)
            
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
                Nica = numel(ica_seq);;
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
                    disp(['before: ',num2str(numel(unique(in_brain_this_ica)))]);
                    
                    in_brain_this_ica = intersect(in_brain,in_brain_this_ica);
                    disp(['  after: ', ...
                          num2str(numel(unique(in_brain_this_ica)))]);
                    
                    all_states(:,j) = mean(base_img(:,in_brain_this_ica),2);
                    
                end
                
                % standard score each dimension separately
                all_states = zscore(all_states);
                
                %% ------------------------------------------------------------
                %% ------------------------------------------------------------
                %% EVC of VALENCE
                
                %% ----------------------------------------
                %% Build tuples for RL
                s_indx = prds.v_indx.h(:,3:(end-1));
                sp_indx = prds.v_indx.h(:,4:end);
                actions = (prds.v_dcmp.h(:,4:end)-prds.v_dcmp.h(:,3:(end-1)));
                errors = (prds.v_dcmp.err(:,4:end));
                terminals = zeros(size(s_indx));terminals(:,end) = 1;
                
                %% reshape to 1D
                s_indx_1d = reshape(s_indx',1,prod(size(s_indx)));
                sp_indx_1d = reshape(sp_indx',1,prod(size(sp_indx)));
                terminals_1d = reshape(terminals',1,prod(size(terminals)));
                actions_1d = zscore(reshape(actions',1,prod(size(actions))));%TICKET
                errors_1d = zscore(reshape(errors',1,prod(size(errors)))); %TICKET
                
                %% assign action to discrete space
                if(act_dscr)
                    %% 5 discrete actions
                    dsc_actions_1d = 0*actions_1d;
                    dsc_actions_1d(find(actions_1d>act_5part(3)))=1; 
                    dsc_actions_1d(find(actions_1d>act_5part(4)))=2; 
                    dsc_actions_1d(find(actions_1d<act_5part(2)))=-1; 
                    dsc_actions_1d(find(actions_1d<act_5part(1)))=-2; 
                else
                    %% 3 discrete actions
                        dsc_actions_1d = 0*actions_1d;
                        dsc_actions_1d(find(actions_1d>act_3part(2)))=1; 
                        dsc_actions_1d(find(actions_1d<act_5part(1)))=-1; 
                end
                
                %% Rescale actions and errors at the group level
                rwd_actions_1d = (dsc_actions_1d-grp_act_v_mean)/grp_act_v_std;
                rwd_errors_1d = (errors_1d-grp_err_v_mean)/grp_err_v_std;
                
                %% Compute reward
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
                cfg.datafile = [proj.path.ctrl.in_evc_opt_mdl,subj_study,'_',name,'_result_v'];
                cfg.maxiter = 1000;
                cfg.regmethod = 'extratrees';
                cfg.singlereg = 0;
                cfg.trees_ntrees = 25;
                cfg.trees_nmin = 2;
                cfg.trees_k = size(X,2)+1;
                cfg.samples = samples;
                
                %% ----------------------------------------
                %% configure fitted Q-iteration
                fittedqi(cfg);
                
            end
            
        end

        %% ----------------------------------------
        %% ----------------------------------------
        %% Analyze control for sign. existence
        %% ----------------------------------------
        %% ----------------------------------------
        sbj_i = 0;

        for i = 1:numel(subjs)
            
            %% extract subject info
            subj_study = subjs{i}.study;
            name = subjs{i}.name;

            try
            
                %% ----------------------------------------
                %% Load subject (state and actions)
                load([proj.path.ctrl.in_evc_opt_mdl,subj_study,'_',name,'_result_v.mat']);
                this_U = Us;
                this_X = Xs';
                
                Q_traj_all = [];
                Q_rand_all = [];
                act_err_all = [];

                cv_ids = setdiff(1:numel(subjs),i);

                for vv = 1:numel(cv_ids)
                    
                    jj = cv_ids(vv);
                    jj

                    %% extract subject info
                    cv_subj_study = subjs{jj}.study;
                    cv_name = subjs{jj}.name;
                    
                    try
                        
                        %% ----------------------------------------
                        %% Load CV subject (Q-functions)
                        load([proj.path.ctrl.in_evc_opt_mdl,cv_subj_study,'_',cv_name,'_result_v.mat']);
                        this_reg = reg;
                        
                        disp(['Loaded CV subject']);
                        
                        Q_traj = [];
                        Q_rand = [];
                        policy = [];
                        opt_policy = [];
                        for k=1:numel(this_U)
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
                        
                        Q_traj_all = [Q_traj_all; Q_traj'];
                        Q_rand_all = [Q_rand_all; Q_rand'];
                        act_err_all = [act_err_all; policy-opt_policy];
                        
                catch
                    disp(' CV subject results not found');
                end
                    
                end %jj
                
                sbj_i = sbj_i+1;
                Q_traj_v_all(a,b,c,sbj_i,:,:) = reshape(median(Q_traj_all),4,30)';
                Q_rand_v_all(a,b,c,sbj_i,:,:) = reshape(median(Q_rand_all),4,30)';
                act_err_v_all(a,b,c,sbj_i,:,:) = reshape(median(act_err_all),4,30)';                

             catch
                 disp('  subject result not found');
             end
            
        end %i
        
        % Save VALENCE intermediate results
        save([proj.path.ctrl.in_evc_opt_mdl,'Q_traj_v_all.mat'],'Q_traj_v_all');
        save([proj.path.ctrl.in_evc_opt_mdl,'Q_rand_v_all.mat'],'Q_rand_v_all');
        save([proj.path.ctrl.in_evc_opt_mdl,'act_err_v_all.mat'],'act_err_v_all');

        % Remove the output of the fitted Q-iteration
        eval(['! rm ',proj.path.ctrl.in_evc_opt_mdl,'*_result_v.mat']);
        
    end
end

toc