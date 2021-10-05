%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

function calc_in_skill(proj,var_name)

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

figure(1)
set(gcf,'color','w');
hold on;

%% ----------------------------------------
%% scatter the underlying stim and feel
predictors = [];
trajs = [];
measures = [];
subjects = [];

sig_cnt = 1;
non_cnt = 1;

b_all = []; %% for power

for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    % log analysis of subject
    logger([subj_study,'_',name],proj.path.logfile);

    try

        %% Load IN trajectory structures
        load([proj.path.ctrl.in_dyn,subj_study,'_',name,'_prds.mat']);

        if(isfield(prds,[var_name,'_dcmp']))
             
            %% extract stims and mean "feel"
            raw_stim = eval(['prds.',var_name,'_dcmp.stim']);
            [~,ids] = sort(raw_stim');
            stim = raw_stim(ids,:);

            feel_raw = eval(['prds.',var_name,'_dcmp.feel(ids,:)']);
            feel_mu = eval(['mean(prds.',var_name,'_dcmp.feel(ids,:),2)']);

            %% remove extreme outliers 
            stim_keep_ids = find(abs(stim)<=3);
            stim_feel_ids = find(abs(feel_mu)<=3);
            cmb_keep_ids = intersect(stim_keep_ids,stim_feel_ids);
            disp(['        # kept stim ids: ', num2str(numel(cmb_keep_ids))]);

            %% extract kept stimuli
            stim_kpt = stim(cmb_keep_ids);
            Nkpt = numel(cmb_keep_ids);
            Nrep = eval(['size(prds.',var_name,'_dcmp.feel,2)']);
            stim_box = repmat(stim_kpt,1,Nrep);
            feel_box = feel_raw(cmb_keep_ids,:);
            traj_box = repmat(1:Nrep,Nkpt,1);

            %% reformat
            stim_form = reshape(stim_box',1,prod(size(stim_box)))';
            feel_form = reshape(feel_box',1,prod(size(feel_box)))';
            traj_form = reshape(traj_box',1,prod(size(traj_box)))';

            stim_form = zscore(stim_form);
            feel_form = zscore(feel_form);

            %% build data for group GLMM
            predictors = [predictors;stim_form]; 
            measures = [measures;feel_form];
            subjects = [subjects;repmat(i,numel(stim_form),1)];
            trajs = [trajs;traj_form];

            %% ----------------------------------------
            %% Quality control below
 
            %% build individual subject structures
            subj = struct();
            subj.study = subj_study;
            subj.name = name;

            tbl = table(feel_form,stim_form,traj_form,'VariableNames', ...
            {'feel','stim','traj'});

            sbj_mdl_fe = fitlme(tbl,['feel ~ 1 + stim + traj']);

            %% Extract Fixed effects
            [~,~,FE] = fixedEffects(sbj_mdl_fe);

            subj.stim = zscore(stim_kpt);
            subj.b1 = FE.Estimate(2); % slope
            subj.b0 = FE.Estimate(1); % intercept
            subj.p1 = FE.pValue(2); %slope
            subj.p0 = FE.pValue(1); %intercept

            %% sort subjects by significance
            if(subj.p1<0.05)
                sig_subjs{sig_cnt} = subj;
                sig_cnt = sig_cnt + 1;
            else
                non_subjs{non_cnt} = subj;
                non_cnt = non_cnt + 1;
            end
            
        else
            disp(['  -Could not find ',var_name,'_dcmp for: ',subj_study,'_', ...
                  name],proj.path.logfile);
        end
        
    catch
        % do nothing
        logger(['  -Could not find/load prds for: ',subj_study,'_', ...
                name],proj.path.logfile);
    end

end

%% ----------------------------------------
%% save out subject groups
save([proj.path.analysis.vr_skill,var_name,'_sig_subjs.mat'],'sig_subjs');
save([proj.path.analysis.vr_skill,var_name,'_non_subjs.mat'],'non_subjs');

%% ----------------------------------------
%% ----------------------------------------
%% Group GLMM fit
measures = double(measures-mean(measures));
predictors = double(predictors-mean(predictors));
trajs = double(trajs-mean(trajs));

tbl = table(measures,predictors,trajs,subjects,'VariableNames', ...
            {'measures','predictors','trajs','subjects'});

mdl_fe = fitlme(tbl,['measures ~ 1 + predictors + trajs']);
mdl_re = fitlme(tbl,['measures ~ 1 + predictors + trajs + (predictors|subjects)']);

fe_vs_re = compare(mdl_fe,mdl_re);

mdl = mdl_fe;
if(fe_vs_re.pValue<0.05);
    mdl=mdl_re;
    logger('  ---Random effects matter',proj.path.logfile);
end

%% Extract Fixed effects
[~,~,FE] = fixedEffects(mdl);

%% Compute effect size
Rsqr = mdl.Rsquared.Adjusted;
Fsqr = Rsqr/(1-Rsqr);
logger(['Rsqr_adj=',num2str(Rsqr)],proj.path.logfile);
logger(['Fsqr=',num2str(Fsqr)],proj.path.logfile);

%% Save out model;
save([proj.path.analysis.vr_skill,var_name,'_vr_skill_mdl.mat'],'mdl');

%% ----------------------------------------
%% ----------------------------------------
%% Plotting

%% scatter raw data
scatter(predictors,measures,10,'MarkerFaceColor', ...
        proj.param.plot.white,'MarkerEdgeColor', ...
        proj.param.plot.light_grey);

%%  ----------------------------------------
%% overlay the individual VR skill plots
for i =1:numel(non_subjs)
    plot(non_subjs{i}.stim,non_subjs{i}.stim*non_subjs{i}.b1+ ...
         non_subjs{i}.b0,'Color',proj.param.plot.light_grey, ...
         'LineWidth',1);
end

for i =1:numel(sig_subjs)
    plot(sig_subjs{i}.stim,sig_subjs{i}.stim*sig_subjs{i}.b1+ ...
         sig_subjs{i}.b0,'Color',proj.param.plot.dark_grey, ...
         'LineWidth',2);
end

%% ----------------------------------------
%% identify max/min x-range|y-rang
%%
xmin = -3;
xmax = 3;
ymin = -2;
ymax = 2;

%% ----------------------------------------
%% overlay VR goal
xseq = linspace(xmin,xmax);
plot(xseq,xseq,'k:','LineWidth',2)

%% ----------------------------------------
%% overlay the group VR skill plot
y_hat = FE.Estimate(1) + FE.Estimate(2)*xseq;
plot(xseq,y_hat,'r-','LineWidth',3);

% No longer plot this b/c position move
% %% ----------------------------------------
% %% indicate goal
% text(1.8,1.7,'\itgoal','FontSize', ...
%      proj.param.plot.axisLabelFontSize-3);

%% ----------------------------------------
%% format figure
xlim([xmin,xmax]);
ylim([ymin,ymax]);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

%% ----------------------------------------
%% explot hi-resolution figure
export_fig skill_summary.png -r300  
eval(['! mv ',proj.path.code,'skill_summary.png ',proj.path.fig,'VR_',var_name,'_skill_summary.png']);

% clean-up
close all;

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

n_tot
n_sig

logger(['----------------------------------------'], ...
       proj.path.logfile);
logger(['Statistical Summary'],proj.path.logfile);
logger(['  -Group effect: ',ge,', p=',gep],proj.path.logfile);
logger(['  -Percent sign. subjs. (p<0.05): ',ssp,'% (',ns,'/',nt, ...
        ')'],proj.path.logfile);
