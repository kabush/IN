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
logger(['Analyzing MVPA GS vs GM                      '],proj.path.logfile);
logger(['*********************************************'],proj.path.logfile);

%% ----------------------------------------
%% load subjs
subjs = proj.process.subjs;

%% ----------------------------------------
%% scatter the underlying stim and feel
v_predictors = [];
v_measures = [];

a_predictors = [];
a_measures = [];

subjects = [];


cnt = 0;
for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    % log analysis of subject
    logger([subj_study,'_',name],proj.path.logfile);

    try

        %% Load IN trajectory structures
        load([proj.path.mvpa.fmri_ex_gs_cls,subj_study,'_',name,'_prds.mat']);
        gs_prds = prds;
        
        load([proj.path.mvpa.fmri_ex_gm_cls,subj_study,'_',name,'_prds.mat']);
        gm_prds = prds;

        v_predictors = [v_predictors;mean(gs_prds.v_cls_hd,1)'];
        v_measures = [v_measures;mean(gm_prds.v_cls_hd,1)'];

        a_predictors = [a_predictors;mean(gs_prds.a_cls_hd,1)'];
        a_measures = [a_measures;mean(gm_prds.a_cls_hd,1)'];
        
        cnt = cnt + 1;
        subjects = [subjects;repmat(cnt,length(gs_prds.v_cls_hd), ...
                                    1)];
        
        
    catch
        % do nothing
        logger(['  -Could not find/load prds for: ',subj_study,'_', ...
                name],proj.path.logfile);
    end

end

%% format figure axes
ymin = -3;
ymax = 3;
xmin = -3;
xmax = 3;


%% ========================================
%% ========================================
%% PLOT VALENCE
%% ========================================
%% ========================================
logger('Analyzing GS vs GM (VALENCE)',proj.path.logfile);

%% format data
v_predictors = double(v_predictors);
v_measures = double(v_measures);
subjects = double(subjects);

%% Group GLMM fit
tbl = table(v_measures,v_predictors,subjects,'VariableNames', ...
            {'trg','pred','subj'});

mdl_fe = fitlme(tbl,['trg ~ 1 + pred']);
mdl_re = fitlme(tbl,['trg ~ 1 + pred + (1+pred|subj)']);

%%Explore random effects across model types
fe_v_re = compare(mdl_fe,mdl_re);

mdl = mdl_fe;
if(fe_v_re.pValue<0.05);
    logger('   random effects matter',proj.path.logfile);
    mdl = mdl_re;
else
    logger('   random effects DO NOT matter',proj.path.logfile);
end

%% Examine Main Effect
[~,~,FE] = fixedEffects(mdl);
if(FE.pValue(2)<0.05)
    logger('   fixed effects are significant',proj.path.logfile);
    logger(['   p=',num2str(FE.pValue(2))],proj.path.logfile);
end

%% compute effect size
Rsqr = mdl.Rsquared.Adjusted;
Fsqr = Rsqr/(1-Rsqr);
logger(['   Rsqr=',num2str(Rsqr)],proj.path.logfile);
logger(['   Fsqr=',num2str(Fsqr)],proj.path.logfile);
logger(' ',proj.path.logfile);

%% plot
figure(1)
set(gcf,'color','w');
scatter(v_predictors,v_measures,10,'MarkerFaceColor', ...
        proj.param.plot.white,'MarkerEdgeColor', ...
        proj.param.plot.very_light_grey);
hold on;


%% overlay the group fixed effect plot
vseq = linspace(xmin,xmax);
y_hat = FE.Estimate(1) + FE.Estimate(2)*vseq;
plot(vseq,y_hat,'r-','LineWidth',3);

xlim([xmin,xmax]);
ylim([ymin,ymax]);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

xlabel('GS-based MVPA Predictions');
ylabel('GM-based MVPA Predictions');

%% ----------------------------------------
%% explot hi-resolution figure
export_fig 'EX_mvpa_gs_vs_gm_val_summary.png' -r300  
eval(['! mv ',proj.path.code,'EX_mvpa_gs_vs_gm_val_summary.png ',proj.path.fig]);

%% ========================================
%% ========================================
%% PLOT AROUSAL
%% ========================================
%% ========================================
logger('Analyzing GS vs GM (AROUSAL)',proj.path.logfile);

%% format data
a_predictors = double(a_predictors);
a_measures = double(a_measures);
subjects = double(subjects);

%% Group GLMM fit
tbl = table(a_measures,a_predictors,subjects,'VariableNames', ...
            {'trg','pred','subj'});

mdl_fe = fitlme(tbl,['trg ~ 1 + pred']);
mdl_re = fitlme(tbl,['trg ~ 1 + pred + (1+pred|subj)']);

%%Explore random effects across model types
fe_a_re = compare(mdl_fe,mdl_re);

mdl = mdl_fe;
if(fe_a_re.pValue<0.05);
    logger('   random effects matter',proj.path.logfile);
    mdl = mdl_re;
else
    logger('   random effects DO NOT matter',proj.path.logfile);
end

%% Examine Main Effect
[~,~,FE] = fixedEffects(mdl);
if(FE.pValue(2)<0.05)
    logger('   fixed effects are significant',proj.path.logfile);
    logger(['   p=',num2str(FE.pValue(2))],proj.path.logfile);
end

%% compute effect size
SS_res=sum((mdl.residuals).^2);
SS_tot=sum((a_measures-mean(a_measures)).^2);
Rsqr = 1-(SS_res/SS_tot);
Fsqr = Rsqr/(1-Rsqr);
logger(['   Rsqr=',num2str(Rsqr)],proj.path.logfile);
logger(['   Fsqr=',num2str(Fsqr)],proj.path.logfile);
logger(' ',proj.path.logfile);

%% plot
figure(2)
set(gcf,'color','w');
scatter(a_predictors,a_measures,10,'MarkerFaceColor', ...
        proj.param.plot.white,'MarkerEdgeColor', ...
        proj.param.plot.very_light_grey);
hold on;


%% ----------------------------------------
%% overlay the group fixed effect plot
vseq = linspace(xmin,xmax);
y_hat = FE.Estimate(1) + FE.Estimate(2)*vseq;
plot(vseq,y_hat,'r-','LineWidth',3);

xlim([xmin,xmax]);
ylim([ymin,ymax]);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

xlabel('GS-based MVPA Predictions');
ylabel('GM-based MVPA Predictions');

%% ----------------------------------------
%% explot hi-resolution figure
export_fig 'EX_mvpa_gs_vs_gm_aro_summary.png' -r300  
eval(['! mv ',proj.path.code,'EX_mvpa_gs_vs_gm_aro_summary.png ',proj.path.fig]);









