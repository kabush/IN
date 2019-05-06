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
logger(['Plotting "Feel" dynamics of Subjects             '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% ----------------------------------------
%% load subjs
load([proj.path.analysis.vr_skill,'a_sig_subjs.mat']);

pos_sbj_id = 21; %fine
neg_sbj_id = 1; %20 not bad
ntr_sbj_id = 2; %fine

%%Statistic on which to measure trajectory
cimed = @(x)median(x);
Nboot = 1000;

max_pos = [];
min_pos = [];
for i=1:size(sig_subjs{pos_sbj_id}.traj,2)
    ci = bootci(Nboot,cimed,sig_subjs{pos_sbj_id}.traj(:,i));
    max_pos = [max_pos,ci(2)];
    min_pos = [min_pos,ci(1)];
end
mu_pos = median(sig_subjs{pos_sbj_id}.traj);

max_neg = [];
min_neg = [];
for i=1:size(sig_subjs{neg_sbj_id}.traj,2)
    ci = bootci(1000,cimed,sig_subjs{neg_sbj_id}.traj(:,i));
    max_neg = [max_neg,ci(2)];
    min_neg = [min_neg,ci(1)];
end
mu_neg = median(sig_subjs{neg_sbj_id}.traj);

max_ntr = [];
min_ntr = [];
for i=1:size(sig_subjs{ntr_sbj_id}.traj,2)
    ci = bootci(1000,cimed,sig_subjs{ntr_sbj_id}.traj(:,i));
    max_ntr = [max_ntr,ci(2)];
    min_ntr = [min_ntr,ci(1)];
end
mu_ntr = median(sig_subjs{ntr_sbj_id}.traj);


%% ----------------------------------------
%% ----------------------------------------
%% ----------------------------------------
%% ----------------------------------------

%% ----------------------------------------
%% overlay VR goal
vseq = linspace(0,size(mu_traj_a,2));
plot(vseq,0*vseq,'k:','LineWidth',2)
hold on;

plot(max_ntr,'Color',proj.param.plot.light_grey,'LineWidth', 1);
plot(min_ntr,'Color',proj.param.plot.light_grey,'LineWidth',1);
plot(mu_ntr,'Color',proj.param.plot.light_grey,'LineWidth',2);

plot(max_pos,'r-','LineWidth',1);
plot(min_pos,'r-','LineWidth',1);
plot(mu_pos,'r-','LineWidth',2);

plot(max_neg,'b-','LineWidth',1);
plot(min_neg,'b-','LineWidth',1);
plot(mu_neg,'b-','LineWidth',2);


%% ----------------------------------------
%% Set y limits
mu_traj_a = [];
for i=1:numel(sig_subjs)
    mu_traj_a = [mu_traj_a;sig_subjs{i}.mu_traj_a];
end
hi_rng = max(max(mu_traj_a))+0.1;
lo_rng = min(min(mu_traj_a))-0.25;

%% ----------------------------------------
%% fill in sections of trajectories
h = fill([0,1,1,0],[lo_rng,lo_rng,hi_rng,hi_rng],'b','edgecolor', ...
         'none');
set(h,'facealpha',.3);
h = fill([1,2,2,1],[lo_rng,lo_rng,hi_rng,hi_rng],'b','edgecolor', ...
         'none');
set(h,'facealpha',.2);
h = fill([6,7,7,6],[lo_rng,lo_rng,hi_rng,hi_rng],'b','edgecolor', ...
         'none');
set(h,'facealpha',.1);

%% ----------------------------------------
%% format figure
xlim([0,size(mu_traj_a,2)]);
ylim([lo_rng,hi_rng]);
hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;
xlabel('Volitional Recall (volume)');
ylabel('Valence(modulate)-Valence(cue)');
text(0.22,lo_rng+0.1,'cue','FontSize', ...
     proj.param.plot.axisLabelFontSize);
text(1.15,lo_rng+0.1,'prep','FontSize', ...
     proj.param.plot.axisLabelFontSize);
text(3.3,lo_rng+0.1,'modulate','FontSize', ...
     proj.param.plot.axisLabelFontSize);
text(6.3,lo_rng+0.1,'ITI','FontSize', ...
     proj.param.plot.axisLabelFontSize);
text(0.22,0.075,'\itgoal','FontSize', ...
     proj.param.plot.axisLabelFontSize-4);


%% ----------------------------------------
%% explot hi-resolution figure
export_fig 'VR_a_dyna_extr_subjs.png' -r300  
eval(['! mv ',proj.path.code,'VR_a_dyna_extr_subjs.png ', ...
      proj.path.fig]);

%% ****************************************
%% TICKET
%% ****************************************
%% Cannot figure out how to use -r<resolution> flag with the
%% functional syntax form of the export_fig command, which is
%% requiring me to to write to local directory and move (above)
