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

%% Storage for analysis
all_Q_traj = [];
all_Q_rand = [];
all_Q_best = []
all_Q_worst = [];

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
        
        % Load ICA masks comprising the state space (and grab activations)
        Nica = 5;
        all_states = zeros(size(base_img,1),Nica);
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

            all_states(:,j) = mean(base_img(:,in_brain_this_ica),2);

        end

        % standard score each dimension separately
        all_states = zscore(all_states);

        %% ----------------------------------------
        %% Build tuples for RL
        s_indx = prds.v_indx.h(:,3:(end-1));
        sp_indx = prds.v_indx.h(:,4:end);
        terminals = zeros(size(s_indx));terminals(:,end) = 1;
        actions = (prds.v_dcmp.h(:,4:end)-prds.v_dcmp.h(:,3:(end-1)));
        errors = (prds.v_dcmp.err(:,4:end));
        
        %% reshape to 1D
        s_indx_1d = reshape(s_indx',1,prod(size(s_indx)));
        sp_indx_1d = reshape(sp_indx',1,prod(size(sp_indx)));
        terminals_1d = reshape(terminals',1,prod(size(terminals)));
        actions_1d = zscore(reshape(actions',1,prod(size(actions))));%TICKET
        errors_1d = zscore(reshape(errors',1,prod(size(errors)))); %TICKET

        % assign a discrete action space %% ***TICKET*** partitions chosen to
        % break actions into evenly split thirds
        dsc_actions_1d = 0*actions_1d;
        dsc_actions_1d(find(actions_1d>.33))=1;        
        dsc_actions_1d(find(actions_1d<-.33))=-1;

        %% ***TICKET*** reward function parameters hardcoded here
        f_cost_act = 1.0;
        f_cost_err = 1.0;
        rewards_1d = -f_cost_act*abs(dsc_actions_1d)-f_cost_err*abs(errors_1d);

        %% ----------------------------------------
        %% format data for writing out structure
        %%   N   // number of tuples
        %%   X   // state
        %%   U   // action
        %%   Xp  // next state
        %%   R   // reward (at Xp)
        %%   T   // Xp terminal state (1) or continuing (0)?
        N = numel(s_indx_1d);
        X = all_states(s_indx_1d,:)';
        Xp = all_states(sp_indx_1d,:)';
        R = rewards_1d;
        T = terminals_1d;
        U = adj_actions_1d;

        % load into struct
        samples = varstostruct('N','X','U','Xp','R','T');

        %% ----------------------------------------
        %% configure fitted Q-iteration
        cfg = struct;
        cfg.run = 1;  % set the run flag;
        cfg.gamma = .98; % discount factor
        cfg.U = unique(U); %[-1 1];
        cfg.datadir = [proj.path.code,'tmp'];
        cfg.datafile = [proj.path.ctrl.in_evc_mdl,subj_study,'_',name,'_result'];
        cfg.maxiter = 1000;
        cfg.regmethod = 'extratrees';
        cfg.singlereg = 1;
        cfg.trees_ntrees = 25;
        cfg.trees_nmin = 2;
        cfg.trees_k = size(X,1)+1;
        cfg.samples = samples;

        %% ----------------------------------------
        %% configure fitted Q-iteration
        fittedqi(cfg);

        %% ----------------------------------------
        %% temporary analysis
        load([cfg.datafile,'.mat']);

        Q_traj = [];
        Q_rand = [];
        Q_best = [];
        Q_worst = [];
        for j=1:N

            Q_traj = [Q_traj;Qp(j,find(cfg.U==Us(j)))];
            Q_rand = [Q_rand;Qp(j,randsample(1:3,1))];

            Qbst = find(Qp(j,:)==max(Qp(j,:)));
            Qwst = find(Qp(j,:)==min(Qp(j,:)));

            if(numel(Qbst)>1)
                Qbst = 1; 
            end

            if(numel(Qwst)>1)
                Qwst = 1; 
            end

            Q_best = [Q_best;Qp(j,Qbst)]; 
            Q_worst = [Q_worst;Qp(j,Qwst)];
        end

        Q_traj = reshape(Q_traj,4,30)';
        Q_rand = reshape(Q_rand,4,30)';
        Q_best = reshape(Q_best,4,30)';
        Q_worst = reshape(Q_worst,4,30)';

        all_Q_traj = [all_Q_traj;mean(Q_traj)];
        all_Q_rand = [all_Q_rand;mean(Q_rand)];
        all_Q_best = [all_Q_best;mean(Q_best)];
        all_Q_worst = [all_Q_worst;mean(Q_worst)];

        %% ----------------------------------------
        %% Local analysis of Qvalues.  Split off
        %% into separate script: ***TICKET***
        figure(1)
        set(gcf,'color','w');
        plot(mean(all_Q_traj),'-r','LineWidth',2);
        hold on;
        plot(mean(all_Q_best),':b');
        plot(mean(all_Q_rand),'-b','LineWidth',2);
        plot(mean(all_Q_worst),':b');
        hold off;
        ylim([-4.5,-1]);
        drawnow;

    end

    export_fig 'VR_Qvalues.png' -r300
    eval(['! mv ',proj.path.code,'VR_Qvalues.png ',proj.path.fig]);
    
end
