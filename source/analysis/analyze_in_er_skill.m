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
logger(['Analyzing ER Skill                               '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.analysis.er_skill]);
    eval(['! rm -rf ',proj.path.analysis.er_skill]);
    disp(['Creating ',proj.path.analysis.er_skill]);
    eval(['! mkdir ',proj.path.analysis.er_skill]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

figure(1)
set(gcf,'color','w');

%% ----------------------------------------
%% scatter the underlying stim and feel
predictors = [];
measures = [];
subjects = [];

sig_cnt = 1;
non_cnt = 1;

for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    % log analysis of subject
    logger([subj_study,'_',name],proj.path.logfile);

    try

        %% Load IN trajectory structures
        load([proj.path.ctrl.in_ctrl,subj_study,'_',name,'_prds.mat']);

        if(isfield(prds,'v_dcmp'))
            
            %% extract stims and mean "feel"
            raw_stim = prds.v_dcmp.stim;
            [~,ids] = sort(raw_stim');
            stim = raw_stim(ids,:);
            feel = mean(prds.v_dcmp.feel(ids,:),2);
            feel_traj = prds.v_dcmp.feel(ids,:);

            %% remove extreme outliers (TICKET hardcoded outliers)
            stim_keep_ids = find(abs(stim)<=3);
            stim_feel_ids = find(abs(feel)<=3);
            cmb_keep_ids = intersect(stim_keep_ids,stim_feel_ids);
            disp(['                  ',num2str(numel(cmb_keep_ids))]);
            stim_clean = stim(cmb_keep_ids);
            feel_clean = feel(cmb_keep_ids);
            feel_traj_clean = feel_traj(cmb_keep_ids);

            %% build data for group GLMM
            predictors = [predictors;double(stim)];
            measures = [measures;double(feel)];
            subjects = [subjects;repmat(i,numel(feel),1)];
            
            %% scatter plot specific points        
            scatter(stim_clean,feel_clean,10,'MarkerFaceColor', ...
                    proj.param.plot.white,'MarkerEdgeColor', ...
                    proj.param.plot.very_light_grey);
            hold on;
            
            %% build individual subject structures
            subj = struct();
            subj.study = subj_study;
            subj.name = name;
            
            subj.stim = stim_clean;
            subj.feel = feel_clean;
            subj.feel_traj = feel_traj_clean;

            [b stat] = robustfit(stim_clean,feel_clean);
            subj.b1 = b(2); % slope
            subj.b0 = b(1); % intercept
            subj.p1 = stat.p(2); %slope
            subj.p0 = stat.p(1); %intercept
 
             %% sort subjects by significance
             if(subj.p1<0.05)
                 sig_subjs{sig_cnt} = subj;
                 sig_cnt = sig_cnt + 1;
             else
                 non_subjs{non_cnt} = subj;
                 non_cnt = non_cnt + 1;
             end
             
        else
            disp(['  -Could not find v_dcmp for: ',subj_study,'_',name],proj.path.logfile);
        end
        
    catch
        % do nothing
        logger(['  -Could not find/load prds for: ',subj_study,'_',name],proj.path.logfile);
    end

end

%% ----------------------------------------
%% save out subject groups
save([proj.path.analysis.er_skill,'sig_subjs.mat'],'sig_subjs');
save([proj.path.analysis.er_skill,'non_subjs.mat'],'non_subjs');

%% ----------------------------------------
%% overlay the individual VR skill plots
for i =1:numel(non_subjs)
    plot(non_subjs{i}.stim,non_subjs{i}.stim*non_subjs{i}.b1+ ...
         non_subjs{i}.b0,'Color',proj.param.plot.light_grey,'LineWidth',1);
end

for i =1:numel(sig_subjs)
    plot(sig_subjs{i}.stim,sig_subjs{i}.stim*sig_subjs{i}.b1+ ...
         sig_subjs{i}.b0,'Color',proj.param.plot.dark_grey,'LineWidth',2);

%%     disp([sig_subjs{i}.study,'_',sig_subjs{i}.name,', b0=', ...
%%           num2str(sig_subjs{i}.b0),' b1=',num2str(sig_subjs{i}.b1)]); 
%%     disp(['          p0=', ...
%%           num2str(sig_subjs{i}.p0),' p1=',num2str(sig_subjs{i}.p1)]); 


end


%% ----------------------------------------
%% identify max/min x-range|y-rang
%%
%% TICKET (automate this decision)
% cmb_subjs = [non_subjs,sig_subjs]
% xrngs = [];
% yrngs = [];
% for i=1:numel(cmb_subjs)
%     xrngs = [xrngs, range(cmb_subjs{i}.stim)];
%     yrngs = [yrngs  range(cmb_subjs{i}.feel)];
% end
% [xmin,xmax] = range(xrngs);
% [ymin,ymax] = range(yrngs);
xmin = -3;
xmax = 3;
ymin = -2;
ymax = 2;

%% ----------------------------------------
%% overlay VR goal
vseq = linspace(xmin,xmax);
plot(vseq,vseq,'k:','LineWidth',2)
hold on;

%% ----------------------------------------
%% Group GLMM fit
tbl = table(measures,predictors,subjects,'VariableNames',{'measures','predictors','subjects'});
mdl = fitlme(tbl,'measures ~ 1 + predictors + (1|subjects) + (predictors-1|subjects)');
[~,~,FE] = fixedEffects(mdl);

%% ----------------------------------------
%% overlay the group VR skill plot
y_hat = FE.Estimate(1) + FE.Estimate(2)*vseq;
plot(vseq,y_hat,'r-','LineWidth',3);

%% ----------------------------------------
%% format figure
xlim([xmin,xmax]);
ylim([ymin,ymax]);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

xlabel('Valence(cue)');
ylabel('Mean Valence(modulate)');

%% ----------------------------------------
%% explot hi-resolution figure
export_fig 'ER_v_skill_summary.png' -r300  
eval(['! mv ',proj.path.code,'ER_v_skill_summary.png ',proj.path.fig]);

%% ****************************************
%% TICKET
%% ****************************************
%% Cannot figure out how to use -r<resolution> flag with the
%% functional syntax form of the export_fig command, which is
%% requiring me to to write to local directory and move (above)

%% ----------------------------------------
%% Output summary

% fixed effects
ge = num2str(FE.Estimate(2));
gep = num2str(FE.pValue(2));

% single subject results
n_sig = numel(sig_subjs);
n_tot = numel(sig_subjs)+numel(non_subjs);
ns = num2str(n_sig);
nt = num2str(n_tot);
ssp = num2str(100*(n_sig/n_tot));

logger(['----------------------------------------'],proj.path.logfile);
logger(['Statistical Summary'],proj.path.logfile);
logger(['  -Group effect: ',ge,', p=',gep],proj.path.logfile);
logger(['  -Percent sign. subjs. (p<0.05): ',ssp,'% (',ns,'/',nt,')'],proj.path.logfile);

