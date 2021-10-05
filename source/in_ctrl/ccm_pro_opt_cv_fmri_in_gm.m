%%========================================
%%========================================
%%
%% Keith Bush, PhD (2019)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

%% Load in path data
load('proj.mat');

%% Initialize log section
logger(['*************************************************'],proj.path.logfile);
logger(['Compute PRO (optional) Trajectories  '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.ctrl.in_pro_opt_mdl]);
    eval(['! rm -rf ',proj.path.ctrl.in_pro_opt_mdl]);
    disp(['Creating ',proj.path.ctrl.in_pro_opt_mdl]);
    eval(['! mkdir ',proj.path.ctrl.in_pro_opt_mdl]);
end


logger(['****************************************'],proj.path.logfile);
logger([' Analyze PRO fits (VALENCE) '],proj.path.logfile);
logger(['****************************************'],proj.path.logfile);

%% ----------------------------------------
%% VALENCE inter-subject PRO predictions
affect_name = 'v';

% Predict the PRO
[pro_all_out,pro_all,mdl] = eval_pro_cv(proj,affect_name);

% Analyze the PRO
[Ntrials,Nsbjs] = size(pro_all_out);
measures = [];
predictors = [];
subjects = [];

sig_cnt = 1;
non_cnt = 1;

for i=1:Nsbjs

    pa = zscore(pro_all(:,i));
    pao = zscore(pro_all_out(:,i));

    measures = [measures;pa];
    predictors = [predictors;pao];%zscore(pro_all_out(:,i))];
    subjects = [subjects;repmat(i,Ntrials,1)];

    % fit individual subjects
    tbl = table(pa,pao,'VariableNames',{'trg','pred'});
    sbj_mdl_fe = fitlme(tbl,['trg ~ 1 + pred']);
    [~,~,FE] = fixedEffects(sbj_mdl_fe);

    % Extract data
    subj = struct();
    subj.stim = pao;
    subj.b1 = FE.Estimate(2);
    subj.b0 = FE.Estimate(1);
    subj.p1 = FE.pValue(2);
    subj.p0 = FE.pValue(1);

    if(subj.p1<0.05)
        sig_subjs{sig_cnt} = subj;
        sig_cnt = sig_cnt + 1;
    else
        non_subjs{non_cnt} = subj;
        non_cnt = non_cnt + 1;
    end

end

xmin = -3;
xmax = 3;
ymin = -2;
ymax = 2;

figure(1)
set(gcf,'color','w');

scatter(predictors,measures,10,'MarkerFaceColor',...
        proj.param.plot.white,'MarkerEdgeColor',...
        proj.param.plot.light_grey);
hold on;

if(non_cnt>1)
    for i=1:numel(non_subjs)
        plot(sort(non_subjs{i}.stim),sort(non_subjs{i}.stim)*non_subjs{i}.b1+non_subjs{i}.b0,'Color',proj.param.plot.light_grey,'LineWidth',1);
    end
end

if(sig_cnt>1)
    for i=1:numel(sig_subjs)
        plot(sort(sig_subjs{i}.stim),sort(sig_subjs{i}.stim)*sig_subjs{i}.b1+sig_subjs{i}.b0,'Color',proj.param.plot.dark_grey,'LineWidth',2);
    end
end

%% Fit main model
tbl = table(measures,predictors,subjects,'VariableNames',{'measures','predictors','subjects'});
mdl_fe = fitlme(tbl,['measures ~ 1 + predictors']);
mdl_re = fitlme(tbl,['measures ~ 1 + predictors + (predictors|subjects)']);
fe_vs_re = compare(mdl_fe,mdl_re);
mdl = mdl_fe;
if(fe_vs_re.pValue<0.05)
    mdl=mdl_re;
    logger('  ---Random effects matters',proj.path.logfile');
end

[~,~,FE] = fixedEffects(mdl);

Rsqr = mdl.Rsquared.Adjusted;
Fsqr = Rsqr/(1-Rsqr);
logger(['Rsqr_adj=',num2str(Rsqr)],proj.path.logfile);
logger(['Fsqr=',num2str(Fsqr)],proj.path.logfile);


xseq = linspace(xmin,xmax);
y_hat = FE.Estimate(1) + FE.Estimate(2)*xseq;
plot(xseq,y_hat,'r-','LineWidth',3);

xlim([xmin,xmax]);
ylim([ymin,ymax]);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

export_fig pro_summary.png -r300
eval(['! mv ',proj.path.code,'pro_summary.png ',proj.path.fig, ...
      'IN_PRO_',affect_name,'_summary.png']);


% clean-up 
close all;

save([proj.path.ctrl.in_pro_opt_mdl,'pro_all_out_',affect_name,'.mat'],'pro_all_out');
save([proj.path.ctrl.in_pro_opt_mdl,'pro_all_',affect_name,'.mat'],'pro_all');
save([proj.path.ctrl.in_pro_opt_mdl,'evaluate_pro_',affect_name,'.mat'],'mdl');



logger(['****************************************'],proj.path.logfile);
logger([' Analyze PRO fits (AROUSAL) '],proj.path.logfile);
logger(['****************************************'],proj.path.logfile);

%% ----------------------------------------
%% AROUSAL inter-subject PRO predictions
affect_name = 'a';

% Predict the PRO
[pro_all_out,pro_all,mdl] = eval_pro_cv(proj,affect_name);

% Analyze the PRO
[Ntrials,Nsbjs] = size(pro_all_out);
measures = [];
predictors = [];
subjects = [];

sig_cnt = 1;
non_cnt = 1;

for i=1:Nsbjs

    pa = zscore(pro_all(:,i));
    pao = zscore(pro_all_out(:,i));

    measures = [measures;pa];
    predictors = [predictors;pao];%zscore(pro_all_out(:,i))];
    subjects = [subjects;repmat(i,Ntrials,1)];

    % fit individual subjects
    tbl = table(pa,pao,'VariableNames',{'trg','pred'});
    sbj_mdl_fe = fitlme(tbl,['trg ~ 1 + pred']);
    [~,~,FE] = fixedEffects(sbj_mdl_fe);

    % Extract data
    subj = struct();
    subj.stim = pao;
    subj.b1 = FE.Estimate(2);
    subj.b0 = FE.Estimate(1);
    subj.p1 = FE.pValue(2);
    subj.p0 = FE.pValue(1);

    if(subj.p1<0.05)
        sig_subjs{sig_cnt} = subj;
        sig_cnt = sig_cnt + 1;
    else
        non_subjs{non_cnt} = subj;
        non_cnt = non_cnt + 1;
    end

end

xmin = -3;
xmax = 3;
ymin = -2;
ymax = 2;

figure(1)
set(gcf,'color','w');

scatter(predictors,measures,10,'MarkerFaceColor',...
        proj.param.plot.white,'MarkerEdgeColor',...
        proj.param.plot.light_grey);
hold on;

if(non_cnt>1)
    for i=1:numel(non_subjs)
        plot(sort(non_subjs{i}.stim),sort(non_subjs{i}.stim)*non_subjs{i}.b1+non_subjs{i}.b0,'Color',proj.param.plot.light_grey,'LineWidth',1);
    end
end

if(sig_cnt>1)
    for i=1:numel(sig_subjs)
        plot(sort(sig_subjs{i}.stim),sort(sig_subjs{i}.stim)*sig_subjs{i}.b1+sig_subjs{i}.b0,'Color',proj.param.plot.dark_grey,'LineWidth',2);
    end
end

%% Fit main model
tbl = table(measures,predictors,subjects,'VariableNames',{'measures','predictors','subjects'});
mdl_fe = fitlme(tbl,['measures ~ 1 + predictors']);
mdl_re = fitlme(tbl,['measures ~ 1 + predictors + (predictors|subjects)']);
fe_vs_re = compare(mdl_fe,mdl_re);
mdl = mdl_fe;
if(fe_vs_re.pValue<0.05)
    mdl=mdl_re;
    logger('  ---Random effects matters',proj.path.logfile');
end

[~,~,FE] = fixedEffects(mdl);

Rsqr = mdl.Rsquared.Adjusted;
Fsqr = Rsqr/(1-Rsqr);
logger(['Rsqr_adj=',num2str(Rsqr)],proj.path.logfile);
logger(['Fsqr=',num2str(Fsqr)],proj.path.logfile);


xseq = linspace(xmin,xmax);
y_hat = FE.Estimate(1) + FE.Estimate(2)*xseq;
plot(xseq,y_hat,'r-','LineWidth',3);

xlim([xmin,xmax]);
ylim([ymin,ymax]);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

export_fig pro_summary.png -r300
eval(['! mv ',proj.path.code,'pro_summary.png ',proj.path.fig, ...
      'IN_PRO_',affect_name,'_summary.png']);

% clean-up 
close all;

save([proj.path.ctrl.in_pro_opt_mdl,'pro_all_out_',affect_name,'.mat'],'pro_all_out');
save([proj.path.ctrl.in_pro_opt_mdl,'pro_all_',affect_name,'.mat'],'pro_all');
save([proj.path.ctrl.in_pro_opt_mdl,'evaluate_pro_',affect_name,'.mat'],'mdl');




