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
load([proj.path.analysis.er_skill,'sig_subjs']);

figure(1)
set(gcf,'color','w');

%% ----------------------------------------
%% scatter the underlying stim and feel
mu_traj_v = [];
sbj_labs = [];

for i = 1:numel(sig_subjs);

    subj_study = sig_subjs{i}.study;
    name = sig_subjs{i}.name;

    % log analysis of subject
    logger([subj_study,'_',name],proj.path.logfile);

    try
        %% Load IN trajectory structures
        load([proj.path.ctrl.in_ctrl,subj_study,'_',name,'_prds.mat']);

        %% Analyze 
        if(isfield(prds,'v_dcmp'))
            
            traj_v = prds.v_dcmp.h-prds.v_dcmp.h(:,1);
            in_ctrl = 1; % set by default to true (search for loss)
            
            %% Determine if subject lost control (TICKET, remove hardcode)
            for j=3:6
                p = signrank(traj_v(:,j));
                if(p<0.05);
                    in_ctrl = 0;
                end
            end

            %% Determine point of max loss of control (TICKET,
            %% remove hardcode))
            if(~in_ctrl)
                p_best = 1.0;
                id_best = -1;
                for j=3:6
                    p = signrank(traj_v(:,j));
                    if(p<p_best)
                        p_best = p;
                        id_best = j;
                    end
                end
            end

            %% Assign subject label
            if(~in_ctrl & p_best<0.05)
                sbj_labs = [sbj_labs; sign(median(traj_v(:,id_best)))];
            else
                sbj_labs = [sbj_labs; 0];
            end
            
            %% Assign median subject trajectory
            mu_traj_v = [mu_traj_v;median(traj_v)];

        else
            disp(['  -Could not find v_dcmp for: ',subj_study,'_',name],proj.path.logfile);
        end

    catch
        % do nothing
        logger(['  -Could not find load prds for: ',subj_study,'_',name],proj.path.logfile);
    end

end

%% ----------------------------------------
%% transfer labels and mu trajectories to
%% structures and save out
for i=1:numel(sig_subjs)
    sig_subjs{i}.in_ctrl = sbj_labs(i);
    sig_subjs{i}.mu_traj_v = mu_traj_v(i,:);
end
save([proj.path.analysis.er_skill,'sig_subjs.mat'],'sig_subjs');

%% ----------------------------------------
%% overlay the individual VR skill plots
ctrl_ids = find(sbj_labs==0);
plot(mu_traj_v(ctrl_ids,:)','Color', proj.param.plot.dark_grey,'LineWidth',2);
hold on;

pos_ids = find(sbj_labs==1);
plot(mu_traj_v(pos_ids,:)','r-','LineWidth',2);
hold on;

neg_ids = find(sbj_labs==-1);
plot(mu_traj_v(neg_ids,:)','b-','LineWidth',2);
hold on;

%% ----------------------------------------
%% overlay VR goal
vseq = linspace(0,size(mu_traj_v,2));
plot(vseq,0*vseq,'k:','LineWidth',2)
hold on;

%% ----------------------------------------
%% Set y limits
hi_rng = 1.0;
lo_rng = -1.0;

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
ylabel('Valence(modulate)-Valence(cue)');
text(0.22,lo_rng+0.1,'cue','FontSize',proj.param.plot.axisLabelFontSize);
text(1.15,lo_rng+0.1,'prep','FontSize',proj.param.plot.axisLabelFontSize);
text(3.3,lo_rng+0.1,'modulate','FontSize',proj.param.plot.axisLabelFontSize);
text(6.3,lo_rng+0.1,'ITI','FontSize',proj.param.plot.axisLabelFontSize);

%% ----------------------------------------
%% explot hi-resolution figure
export_fig 'ER_v_dyna_sbj_labels.png' -r300  
eval(['! mv ',proj.path.code,'ER_v_dyna_sbj_labels.png ',proj.path.fig]);

%% ****************************************
%% TICKET
%% ****************************************
%% Cannot figure out how to use -r<resolution> flag with the
%% functional syntax form of the export_fig command, which is
%% requiring me to to write to local directory and move (above)
