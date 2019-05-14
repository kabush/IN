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
logger(['*********************************************'],proj.path.logfile);
logger(['Analyzing EMG responses to IN stimuli        '],proj.path.logfile);
logger(['*********************************************'],proj.path.logfile);

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% ----------------------------------------
%% scatter the underlying stim and feel
zygo_measures = [];
corr_measures = [];
zygo_predictors = [];
corr_predictors = [];
subjects = [];

for i = 1:22 % numel(subjs) ONLY CTM subjs have EMG

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    % debug
    logger([subj_study,'_',name],proj.path.logfile);

    try
        load([proj.path.betas.emg_in_beta,subj_study,'_',name,'_in_betas.mat']);
        load([proj.path.betas.emg_in_beta,subj_study,'_',name,'_feel_betas.mat']);
    catch
        logger('    Could not find emg beta file for processing.',proj.path.logfile);
    end

    in_zygo_betas = zscore([in_betas.zygo.id1;in_betas.zygo.id2]);
    in_corr_betas = zscore([in_betas.corr.id1;in_betas.corr.id2]);

    feel_zygo_betas = [feel_betas.zygo.id1,feel_betas.zygo.id2];
    feel_corr_betas = [feel_betas.corr.id1,feel_betas.corr.id2];

    Nfact= numel(feel_zygo_betas)/numel(in_zygo_betas);
    mu_feel_zygo_betas = zscore(mean(reshape(feel_zygo_betas,Nfact,numel(in_zygo_betas))',2));

    Nfact= numel(feel_corr_betas)/numel(in_corr_betas);
    mu_feel_corr_betas = zscore(mean(reshape(feel_corr_betas,Nfact,numel(in_corr_betas))',2));

    zygo_measures = [zygo_measures;mu_feel_zygo_betas];
    corr_measures = [corr_measures;mu_feel_corr_betas];    

    zygo_predictors = [zygo_predictors,in_zygo_betas'];
    corr_predictors = [corr_predictors,in_corr_betas'];

    subjects = [subjects;repmat(i,numel(in_zygo_betas),1)];

end

disp('----------------------------------------');
disp(' ZYGOMATICUS ANALYSIS ');
disp('----------------------------------------');
%% ----------------------------------------
%% Group GLMM fit
measures = double(zygo_measures);
predictors = double(zygo_predictors');
subjects = double(subjects);

tbl = table(measures,predictors,subjects,'VariableNames',{'trg', ...
                    'pred','subj'});

mdl_fe = fitlme(tbl,['trg ~ 1 + pred']);
mdl_re= fitlme(tbl,['trg ~ 1 + pred + (1+pred|subj)']);

%%Explore random effects across model types
fe_v_re = compare(mdl_fe,mdl_re);

mdl = mdl_fe;
disp('Testing Random Effects');
if(fe_v_re.pValue<0.05);
    disp('  random effects matter');
    mdl = mdl_re;
else
    disp('  random effects DO NOT matter');
end


disp(' ');

%% ----------------------------------------
%% Examine Main Effect
[~,~,FE] = fixedEffects(mdl);
if(FE.pValue(2)<0.05)
    disp('Fixed Effects are significant');
    disp(['  p=',num2str(FE.pValue(2))]);
end


%% ----------------------------------------
%% compute effect size
SS_res=sum((mdl.residuals).^2);
SS_tot=sum((measures-mean(measures)).^2);
Rsqr = 1-(SS_res/SS_tot);
Fsqr = Rsqr/(1-Rsqr);
logger(['  Rsqr=',num2str(Rsqr)],proj.path.logfile);
logger(['  Fsqr=',num2str(Fsqr)],proj.path.logfile);

disp(' ');

figure(1)
set(gcf,'color','w');

%% ----------------------------------------
%% plot all the datapoints
scatter(predictors,measures,60,'MarkerFaceColor', ...
        proj.param.plot.white,'MarkerEdgeColor', ...
        proj.param.plot.dark_grey);
hold on;

%% ----------------------------------------
%% format figure
ymin = -2; 
ymax = 2;
xmin = -2; 
xmax = 2;

%% ----------------------------------------
%% overlay the group VR skill plot
vseq = linspace(xmin,xmax);
y_hat = FE.Estimate(1) + FE.Estimate(2)*vseq;
plot(vseq,y_hat,'k-','LineWidth',6);

xlim([xmin,xmax]);
ylim([ymin,ymax]);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

% xlabel('VR(cue) EMG responses');
% ylabel('VR(modulate) EMG responses');

%% ----------------------------------------
%% explot hi-resolution figure
export_fig 'IN_emg_zygo_summary.png' -r300  
eval(['! mv ',proj.path.code,'IN_emg_zygo_summary.png ',proj.path.fig]);


disp('----------------------------------------');
disp(' CORRUGATOR ANALYSIS ');
disp('----------------------------------------');
%% ----------------------------------------
%% Group GLMM fit
measures = double(corr_measures);
predictors = double(corr_predictors');
subjects = double(subjects);

tbl = table(measures,predictors,subjects,'VariableNames',{'trg', ...
                    'pred','subj'});

mdl_fe = fitlme(tbl,['trg ~ 1 + pred']);
mdl_re= fitlme(tbl,['trg ~ 1 + pred + (1+pred|subj)']);

%%Explore random effects across model types
fe_v_re = compare(mdl_fe,mdl_re);

mdl = mdl_fe;
disp('Testing Random Effects');
if(fe_v_re.pValue<0.05);
    disp('  random effects matter');
    mdl = mdl_re;
else
    disp('  random effects DO NOT matter');
end


disp(' ');

%% ----------------------------------------
%% Examine Main Effect
[~,~,FE] = fixedEffects(mdl);
if(FE.pValue(2)<0.05)
    disp('Fixed Effects are significant');
    disp(['  p=',num2str(FE.pValue(2))]);
end


%% ----------------------------------------
%% compute effect size
SS_res=sum((mdl.residuals).^2);
SS_tot=sum((measures-mean(measures)).^2);
Rsqr = 1-(SS_res/SS_tot);
Fsqr = Rsqr/(1-Rsqr);
logger(['  Rsqr=',num2str(Rsqr)],proj.path.logfile);
logger(['  Fsqr=',num2str(Fsqr)],proj.path.logfile);

disp(' ');

figure(1)
set(gcf,'color','w');

%% ----------------------------------------
%% plot all the datapoints
scatter(predictors,measures,60,'MarkerFaceColor', ...
        proj.param.plot.white,'MarkerEdgeColor', ...
        proj.param.plot.dark_grey);
hold on;

%% ----------------------------------------
%% format figure
ymin = -2; 
ymax = 2;
xmin = -2; 
xmax = 2;

%% ----------------------------------------
%% overlay the group VR skill plot
vseq = linspace(xmin,xmax);
y_hat = FE.Estimate(1) + FE.Estimate(2)*vseq;
plot(vseq,y_hat,'k-','LineWidth',6);

xlim([xmin,xmax]);
ylim([ymin,ymax]);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

% xlabel('VR(cue) EMG responses');
% ylabel('VR(modulate) EMG responses');

%% ----------------------------------------
%% explot hi-resolution figure
export_fig 'IN_emg_corr_summary.png' -r300  
eval(['! mv ',proj.path.code,'IN_emg_corr_summary.png ',proj.path.fig]);

%% ***TICKET DO BOTH ZYGO AND CORR
