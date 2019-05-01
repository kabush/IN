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
logger(['*************************************************'], ...
       proj.path.logfile);
logger(['Plotting "Feel" dynamics                         '], ...
       proj.path.logfile);
logger(['*************************************************'], ...
       proj.path.logfile);

%% ----------------------------------------
%% load subjs
load([proj.path.analysis.vr_skill,'a_sig_subjs.mat']);

figure(1)
set(gcf,'color','w');

%% ----------------------------------------
%% scatter the underlying stim and feel
mu_traj_a = [];
sbj_labs = [];

for i = 1:numel(sig_subjs);

    subj_study = sig_subjs{i}.study;
    name = sig_subjs{i}.name;

    % log analysis of subject
    logger([subj_study,'_',name],proj.path.logfile);


    try
        %% Load IN trajectory structures
        load([proj.path.ctrl.in_ctrl,subj_study,'_',name, ...
              '_prds.mat']);

        %% Analyze 
        if(isfield(prds,'a_dcmp'))
            
            traj_a = prds.a_dcmp.h-prds.a_dcmp.h(:,1);

            % conduct interval analysis
            cimed = @(x)median(x);
            Nboot = 1000;
            max_traj = [];
            min_traj = [];
            for j=1:size(traj_a,2)
                ci = bootci(Nboot,cimed,traj_a(:,j));
                max_traj = [max_traj,ci(2)];
                min_traj = [min_traj,ci(1)];
            end

            neg_ctrl_err = 0-max_traj(3:6);
            pos_ctrl_err = min_traj(3:6);

            max_neg_ctrl_err = neg_ctrl_err(find(neg_ctrl_err==max(neg_ctrl_err)));
            max_pos_ctrl_err = pos_ctrl_err(find(pos_ctrl_err==max(pos_ctrl_err)));

            in_ctrl = 1;
            if(max_neg_ctrl_err <0 & max_pos_ctrl_err <0)
                %% in control over range
                %% do nothing
                sbj_labs = [sbj_labs,0];
            else

                if(max_neg_ctrl_err > max_pos_ctrl_err)
                    %% negative error label
                    sbj_labs = [sbj_labs,-1];
                else
                    %% positive error label
                    sbj_labs = [sbj_labs,1];
                end
            end
%                 

%             in_ctrl = 1; % set by default to true (search for loss)
%             
%             %% Determine if subject lost control (TICKET, remove
%             %% hardcode)
%             for j=3:6
%                 
%                 p = signrank(traj_a(:,j));
%                 if(p<0.05);
%                     in_ctrl = 0;
%                 end
%             end
% 
%             %% Determine point of max loss of control (TICKET,
%             %% remove hardcode))
%             if(~in_ctrl)
%                 p_best = 1.0;
%                 id_best = -1;
%                 for j=3:6
%                     p = signrank(traj_a(:,j));
%                     if(p<p_best)
%                         p_best = p;
%                         id_best = j;
%                     end
%                 end
%             end
% 
%             %% Assign subject label
%             if(~in_ctrl & p_best<0.05)
%                 sbj_labs = [sbj_labs; sign(median(traj_a(:, ...
%                                                          id_best)))];
%             else
%                 sbj_labs = [sbj_labs; 0];
%             end
            
            %% Assign median subject trajectory for visualization
            mu_traj_a = [mu_traj_a;median(traj_a)];

        else
            disp(['  -Could not find a_dcmp for: ',subj_study,'_', ...
                  name],proj.path.logfile);
        end

    catch
        % do nothing
        logger(['  -Could not find load prds for: ',subj_study,'_', ...
                name],proj.path.logfile);
    end

end

%% ----------------------------------------
%% Report out breakdown
logger(['  Control: ',num2str(numel(find(sbj_labs==0))),'/',num2str(numel(sig_subjs))]);
logger(['  Pos out: ',num2str(numel(find(sbj_labs==1))),'/',num2str(numel(sig_subjs))]);
logger(['  Neg out: ',num2str(numel(find(sbj_labs==-1))),'/',num2str(numel(sig_subjs))]);

%% ----------------------------------------
%% transfer labels and mu trajectories to
%% structures and save out
for i=1:numel(sig_subjs)
    sig_subjs{i}.in_ctrl = sbj_labs(i);
    sig_subjs{i}.mu_traj_a = mu_traj_a(i,:);
end
save([proj.path.analysis.vr_skill,'a_sig_subjs.mat'],'sig_subjs');

%% ----------------------------------------
%% overlay the individual VR skill plots
ctrl_ids = find(sbj_labs==0);
plot(mu_traj_a(ctrl_ids,:)','Color', proj.param.plot.dark_grey, ...
     'LineWidth',2);
hold on;

pos_ids = find(sbj_labs==1);
plot(mu_traj_a(pos_ids,:)','r-','LineWidth',2);
hold on;

neg_ids = find(sbj_labs==-1);
plot(mu_traj_a(neg_ids,:)','b-','LineWidth',2);
hold on;

%% ----------------------------------------
%% overlay VR goal
vseq = linspace(0,size(mu_traj_a,2));
plot(vseq,0*vseq,'k:','LineWidth',2)
hold on;

%% ----------------------------------------
%% Set y limits
hi_rng = max(max(mu_traj_a))+0.1;
lo_rng = min(min(mu_traj_a))-0.25;

% %% ----------------------------------------
% %% fill in sections of trajectories
% h = fill([0,1,1,0],[lo_rng,lo_rng,hi_rng,hi_rng],'b','edgecolor', ...
%          'none');
% set(h,'facealpha',.3);
% h = fill([1,2,2,1],[lo_rng,lo_rng,hi_rng,hi_rng],'b','edgecolor', ...
%          'none');
% set(h,'facealpha',.2);
% h = fill([6,7,7,6],[lo_rng,lo_rng,hi_rng,hi_rng],'b','edgecolor', ...
%          'none');
% set(h,'facealpha',.1);

%% ----------------------------------------
%% fill in sections of trajectories ***TICKET*** changed color
%% for CNS poster
cspec = [0.15,0.45,1.0];
h = fill([0,1,1,0],[lo_rng,lo_rng,hi_rng,hi_rng],cspec,'edgecolor', ...
         'none');
set(h,'facealpha',.3);
h = fill([1,2,2,1],[lo_rng,lo_rng,hi_rng,hi_rng],cspec,'edgecolor', ...
         'none');
set(h,'facealpha',.2);
h = fill([6,7,7,6],[lo_rng,lo_rng,hi_rng,hi_rng],cspec,'edgecolor', ...
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
ylabel('Arousal(modulate)-Arousal(cue)');
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
export_fig 'VR_a_dyna_sbj_labels.png' -r300  
eval(['! mv ',proj.path.code,'VR_a_dyna_sbj_labels.png ', ...
      proj.path.fig]);

%% ****************************************
%% TICKET
%% ****************************************
%% Cannot figure out how to use -r<resolution> flag with the
%% functional syntax form of the export_fig command, which is
%% requiring me to to write to local directory and move (above)
