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
logger(['Analyzing VR Q values                             '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% Storage for analysis
% valence
all_v_Q_traj = [];
all_v_Q_rand = [];
all_v_Q_best = []
all_v_Q_worst = [];

% arousal
all_a_Q_traj = [];
all_a_Q_rand = [];
all_a_Q_best = []
all_a_Q_worst = [];

for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    % log analysis of subject
    logger([subj_study,'_',name],proj.path.logfile);

    try

        %% ----------------------------------------
        %% VALENCE analysis
        load([proj.path.ctrl.in_evc_mdl,subj_study,'_',name,'_result_v.mat']);

        Q_traj = [];
        Q_rand = [];
        Q_best = [];
        Q_worst = [];
        for j=1:size(Xs,2)

            Q_traj = [Q_traj;Qp(j,find(cfg.U==Us(j)))];


            Q_lo = Qp(j,1);
            Q_mi = Qp(j,2);
            Q_hi = Qp(j,3);


            Q_rand = [Q_rand;mean([Q_lo,Q_mi,Q_hi])]; %Qp(j,randsample(1:3,1))];

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

        all_v_Q_traj = [all_v_Q_traj;mean(Q_traj)];
        all_v_Q_rand = [all_v_Q_rand;mean(Q_rand)];
        all_v_Q_best = [all_v_Q_best;mean(Q_best)];
        all_v_Q_worst = [all_v_Q_worst;mean(Q_worst)];

        %% ----------------------------------------
        %% AROUSAL analysis
        load([proj.path.ctrl.in_evc_mdl,subj_study,'_',name,'_result_a.mat']);
        
        Q_traj = [];
        Q_rand = [];
        Q_best = [];
        Q_worst = [];
        for j=1:size(Xs,2)
        
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
        
        all_a_Q_traj = [all_a_Q_traj;mean(Q_traj)];
        all_a_Q_rand = [all_a_Q_rand;mean(Q_rand)];
        all_a_Q_best = [all_a_Q_best;mean(Q_best)];
        all_a_Q_worst = [all_a_Q_worst;mean(Q_worst)];

    catch
        disp('   file not found');

    end

end

%% ----------------------------------------
%% ----------------------------------------
%% Plot Group Q-values via time (VALENCE)
figure(1)
set(gcf,'color','w');

%plot
plot(mean(all_v_Q_traj),'-r','LineWidth',2);
hold on;
plot(mean(all_v_Q_best),':b');
plot(mean(all_v_Q_rand),'-b','LineWidth',2);
plot(mean(all_v_Q_worst),':b');
hold off;

% scale axes
ymax = max(mean(all_v_Q_best))+0.1;
ymin = min(mean(all_v_Q_worst))-0.1;
ylim([ymin,ymax]);

%export
export_fig 'VR_Qvalues_v.png' -r300
eval(['! mv ',proj.path.code,'VR_Qvalues_v.png ',proj.path.fig]);

%% ----------------------------------------
%% Statistical tests
p1 = ranksum(all_v_Q_traj(:,1),all_v_Q_rand(:,1));
p2 = ranksum(all_v_Q_traj(:,2),all_v_Q_rand(:,2));
p3 = ranksum(all_v_Q_traj(:,3),all_v_Q_rand(:,3));
p4 = ranksum(all_v_Q_traj(:,4),all_v_Q_rand(:,4));

logger(['----------------------------------------'],proj.path.logfile);
logger(['VALENCE                          '],proj.path.logfile);
logger(['Statistical Summary: traj != rand'],proj.path.logfile);
logger(['   p(t1)=',num2str(p1)],proj.path.logfile);
logger(['   p(t2)=',num2str(p2)],proj.path.logfile);
logger(['   p(t3)=',num2str(p3)],proj.path.logfile);
logger(['   p(t4)=',num2str(p4)],proj.path.logfile);
logger([' ']);


%% ----------------------------------------
%% ----------------------------------------
%% Plot Group Q-values via time (AROUSAL)
figure(2)
set(gcf,'color','w');

%plot
plot(mean(all_a_Q_traj),'-r','LineWidth',2);
hold on;
plot(mean(all_a_Q_best),':b');
plot(mean(all_a_Q_rand),'-b','LineWidth',2);
plot(mean(all_a_Q_worst),':b');
hold off;

% scale axes
ymax = max(mean(all_a_Q_best))+0.1;
ymin = min(mean(all_a_Q_worst))-0.1;
ylim([ymin,ymax]);

%export
export_fig 'VR_Qvalues_a.png' -r300
eval(['! mv ',proj.path.code,'VR_Qvalues_a.png ',proj.path.fig]);

%% ----------------------------------------
%% Statistical tests
p1 = ranksum(all_a_Q_traj(:,1),all_a_Q_rand(:,1));
p2 = ranksum(all_a_Q_traj(:,2),all_a_Q_rand(:,2));
p3 = ranksum(all_a_Q_traj(:,3),all_a_Q_rand(:,3));
p4 = ranksum(all_a_Q_traj(:,4),all_a_Q_rand(:,4));

logger(['----------------------------------------'],proj.path.logfile);
logger(['AROUSAL                          '],proj.path.logfile);
logger(['Statistical Summary: traj != rand'],proj.path.logfile);
logger(['   p(t1)=',num2str(p1)],proj.path.logfile);
logger(['   p(t2)=',num2str(p2)],proj.path.logfile);
logger(['   p(t3)=',num2str(p3)],proj.path.logfile);
logger(['   p(t4)=',num2str(p4)],proj.path.logfile);
logger([' ']);


