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
logger(['Analyzing VR Q values (GRID SEARCH)              '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% Meta RL Parameter search
action_5dscr_set = [1,0]; % otherwise 3 discrete actions
discount_set = [0:.1:1];
reward_frac_set = [0:.2:1]; % balance between reward/action

% Save VALENCE intermediate results
load([proj.path.ctrl.in_evc_opt_mdl,'Q_traj_v_all.mat']);
load([proj.path.ctrl.in_evc_opt_mdl,'Q_rand_v_all.mat']);
load([proj.path.ctrl.in_evc_opt_mdl,'act_err_v_all.mat']);

%% Calculate set-sizes
[Nact,Ndsct,Nfrac,Nsbj,Nvr,Nt] = size(Q_traj_v_all);

%% META-LOOPS HERE
all_q_traj_dsct = zeros(Ndsct,Nfrac,Nsbj);
all_q_rand_dsct = zeros(Ndsct,Nfrac,Nsbj);
all_q_diff_dsct = zeros(Ndsct,Nfrac,Nsbj);
all_q_perf_dsct = zeros(Ndsct,Nfrac,Nsbj);
act_err_dsct = zeros(Ndsct,Nfrac,Nsbj);

for a=1:1 %Nact

    for b=1:Ndsct

        for c=1:Nfrac
            
            act_dscr = action_5dscr_set(a);
            gamma = discount_set(b);
            rwrd_f = reward_frac_set(c);

            all_q_diff = [];
            all_q_rand = [];
            all_q_effx = [];

            for d = 1:Nsbj

                %% Q analysis
                q_traj = squeeze(Q_traj_v_all(a,b,c,d,:,:));
                q_rand = squeeze(Q_rand_v_all(a,b,c,d,:,:));

                q_diff = (q_traj-q_rand);
                q_perf = q_diff./abs(q_rand);
                q_effx = q_diff/std(vec(q_diff));

                all_q_traj_dsct(b,c,d) = median(median(q_traj));
                all_q_rand_dsct(b,c,d) = median(median(q_rand));
                all_q_diff_dsct(b,c,d) = median(median(q_diff));
                all_q_perf_dsct(b,c,d) = median(median(q_perf));

                %% Action analysis
                act_err_dsct(b,c,d) = mean(mean(abs(squeeze(act_err_v_all(a,b,c,d,:,:)))));

            end

            disp(['act=',num2str(act_dscr),', gamma=',num2str(gamma),', frac=',num2str(rwrd_f)]);            

        end
    end
end

%% Observe best parameter combination
all_q_perf = zeros(Ndsct,Nfrac);
all_act_err = zeros(Ndsct,Nfrac);
all_sig_test = zeros(Ndsct,Nfrac);
for b=1:Ndsct
    for c=1:Nfrac
        all_q_perf(b,c) = median(squeeze(all_q_perf_dsct(b,c,:)));
        all_act_err(b,c) = mean(squeeze(act_err_dsct(b,c,:)));
        all_sig_test(b,c) = signrank(squeeze(all_q_perf(b,c,:)));
    end
end

% 
% %% ----------------------------------------
% %% Analyze Q-values (on-policy vs random)
% %% ----------------------------------------
% figure(1)
% set(gcf,'color','w');
% q_perf_plot = []
% for i = 1:Ndsct
%     q_perf_plot = [q_perf_plot,median(median(squeeze(all_q_perf_dsct(i,:,:))))];
% end
% plot(q_perf_plot,'b-','LineWidth',3);
% xticklabels({'0',' ','.2',' ','.4',' ','.6',' ','.8',' ','1'});
% box off;
% hold off;
% fig = gcf;
% ax = fig.CurrentAxes;
% ax.FontSize = proj.param.plot.axisLabelFontSize;
% 
% export_fig 'Q_value_sig_IN_task.png' -r300
% eval(['! mv ',proj.path.code,'Q_value_sig_IN_task.png ',proj.path.fig]);

% %% --------------------------------------------
% %% Analyze Action Errors (on-policy vs optimal)
% %% --------------------------------------------
% figure(2);
% set(gcf,'color','w');
% load([proj.path.ctrl.in_evc_mdl,'act_err_v_all.mat']);
% err = [];
% for i=1:Ndsct
%     err = [err,median(median(squeeze(act_err(i,:))))];
% end
% plot(1:11,err,'b-','LineWidth',3);
% xticklabels({'0',' ','.2',' ','.4',' ','.6',' ','.8',' ','1'});
% box off;
% hold off;
% fig = gcf;
% ax = fig.CurrentAxes;
% ax.FontSize = proj.param.plot.axisLabelFontSize;
% 
% %% Significance test of Q-values
% sig_test = zeros(Ndsct,1);
% for i=1:Ndsct
%     data = squeeze(all_q_perf_dsct(i,:,:));
%     q_sig_test(i) = signrank(reshape(data,[],1));
% end
% 
% %% Significance test of Actions
% sig_test = zeros(Ndsct,1);
% for i=1:Ndsct
%     data1 = squeeze(act_err(2,:));
%     data2 = squeeze(act_err(i,:));
%     act_sig_test(i) = ranksum(data1,data2);
% end
% 
% export_fig 'action_err_IN_task.png' -r300
% eval(['! mv ',proj.path.code,'action_err_IN_task.png ',proj.path.fig]);
