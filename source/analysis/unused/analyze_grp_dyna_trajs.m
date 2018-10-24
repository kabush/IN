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
logger(['Plotting "Feel" dynamics                         '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

figure(1)
set(gcf,'color','w');

%% ----------------------------------------
%% scatter the underlying stim and feel
mu_traj_v = [];
max_traj_dv = [];

for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    % log analysis of subject
    logger([subj_study,'_',name],proj.path.logfile);

    try
        %% Load IN trajectory structures
        %load([proj.path.ctrl.in_ctrl_loocv,subj_study,'_',name,'_prds_loocv.mat']);
        load([proj.path.ctrl.in_ctrl,subj_study,'_',name,'_prds.mat']);
    catch
        % do nothing
        logger(['  -Could not find load prds for: ',subj_study,'_',name],proj.path.logfile);
    end

    if(isfield(prds,'v_dcmp'))

        traj_v = prds.v_dcmp.h-prds.v_dcmp.h(:,1);
        mu_traj_v = [mu_traj_v;median(traj_v)];

    else
        disp(['  -Could not find v_dcmp for: ',subj_study,'_',name],proj.path.logfile);
    end

end

%% ----------------------------------------
%% overlay the individual VR skill plots
for i = 1:size(mu_traj_v,1)
    plot(1:size(mu_traj_v,2),mu_traj_v(i,:),'Color',proj.param.plot.light_grey,'LineWidth',2);
    hold on;
end

%% ----------------------------------------
%% overlay the mean CI trajectory
hi_vec = [];
lo_vec = [];
for i = 1:size(mu_traj_v,2)
    [h p ci stat] = ttest(mu_traj_v(:,i));
    hi_vec = [hi_vec,ci(2)];
    lo_vec = [lo_vec,ci(1)];
end

plot(1:size(mu_traj_v,2),hi_vec,'r-','LineWidth',3);
plot(1:size(mu_traj_v,2),lo_vec,'r-','LineWidth',3);

%% ----------------------------------------
%% overlay VR goal
vseq = linspace(0,size(mu_traj_v,2));
plot(vseq,0*vseq,'k:','LineWidth',2)
hold on;

%% ----------------------------------------
%% Set y limits
hi_rng = 1.25;
lo_rng = -1.5;

%% ----------------------------------------
%% fill in sections of trajectories
h = fill([0,1,1,0],[lo_rng,lo_rng,hi_rng,hi_rng],'b','edgecolor','none');
set(h,'facealpha',.3);
h = fill([1,2,2,1],[lo_rng,lo_rng,hi_rng,hi_rng],'b','edgecolor','none');
set(h,'facealpha',.2);
h = fill([6,7,7,6],[lo_rng,lo_rng,hi_rng,hi_rng],'b','edgecolor','none');
set(h,'facealpha',.1);

%% ----------------------------------------
%% format figure
xlim([0,size(mu_traj_v,2)]);
ylim([lo_rng,hi_rng]);
hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;
xlabel('Volitional Recall (volume)');
ylabel('Valence(mod)-Valence(cue)');
text(0.22,lo_rng+0.1,'cue','FontSize',proj.param.plot.axisLabelFontSize);
text(1.15,lo_rng+0.1,'prep','FontSize',proj.param.plot.axisLabelFontSize);
text(3.3,lo_rng+0.1,'modulate','FontSize',proj.param.plot.axisLabelFontSize);
text(6.3,lo_rng+0.1,'ITI','FontSize',proj.param.plot.axisLabelFontSize);

%% ----------------------------------------
%% explot hi-resolution figure
export_fig 'ER_v_dyna_grp_summary.png' -r300  
eval(['! mv ',proj.path.code,'ER_v_dyna_summary.png ',proj.path.fig]);

%% ****************************************
%% TICKET
%% ****************************************
%% Cannot figure out how to use -r<resolution> flag with the
%% functional syntax form of the export_fig command, which is
%% requiring me to to write to local directory and move (above)
