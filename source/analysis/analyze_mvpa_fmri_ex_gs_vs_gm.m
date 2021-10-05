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

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.analysis.ex_gs_vs_gm]);
    eval(['! rm -rf ',proj.path.analysis.ex_gs_vs_gm]);
    disp(['Creating ',proj.path.analysis.ex_gs_vs_gm]);
    eval(['! mkdir ',proj.path.analysis.ex_gs_vs_gm]);
end

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
v_sig_cnt = 1;
v_non_cnt = 1;

a_predictors = [];
a_measures = [];
a_sig_cnt = 1;
a_non_cnt = 1;

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

        %% Load EX predictionsry structures
        load([proj.path.mvpa.fmri_ex_gs_cls,subj_study,'_',name,'_prds.mat']);
        gs_prds = prds;
        
        load([proj.path.mvpa.fmri_ex_gm_cls,subj_study,'_',name,'_prds.mat']);
        gm_prds = prds;

        %% ----------------------------------------
        v_predictors = [v_predictors;zscore(mean(gs_prds.v_cls_hd,1))'];
        v_measures = [v_measures;zscore(mean(gm_prds.v_cls_hd,1))'];

        % model subject level
        subj = struct();
        subj.study = subj_study;
        subj.name = name;
        
        pred = double(zscore(mean(gs_prds.v_cls_hd,1))');
        trg = double(zscore(mean(gm_prds.v_cls_hd,1))');
        tbl = table(pred,trg,'VariableNames',{'pred','trg'});

        sbj_mdl_fe = fitlme(tbl,['trg ~ 1 + pred']);
 
        % extract fixed effects
        [~,~,FE] = fixedEffects(sbj_mdl_fe);

        subj.stim = pred;
        subj.b1 = FE.Estimate(2);
        subj.b0 = FE.Estimate(1);
        subj.p1 = FE.pValue(2);
        subj.p0 = FE.pValue(1);

        if(subj.p1<0.05)
            v_sig_subjs{v_sig_cnt} = subj;
            v_sig_cnt = v_sig_cnt + 1;
        else
            v_non_subjs{v_non_cnt} = subj;
            v_non_cnt = v_non_cnt +1;
        end


        %% ----------------------------------------
        a_predictors = [a_predictors;zscore(mean(gs_prds.a_cls_hd,1))'];
        a_measures = [a_measures;zscore(mean(gm_prds.a_cls_hd,1))'];

        % model subject level
        subj = struct();
        subj.study = subj_study;
        subj.name = name;
        
        pred = double(zscore(mean(gs_prds.a_cls_hd,1))');
        trg = double(zscore(mean(gm_prds.a_cls_hd,1))');
        tbl = table(pred,trg,'VariableNames',{'pred','trg'});

        sbj_mdl_fe = fitlme(tbl,['trg ~ 1 + pred']);
 
        % extract fixed effects
        [~,~,FE] = fixedEffects(sbj_mdl_fe);

        subj.stim = pred;
        subj.b1 = FE.Estimate(2);
        subj.b0 = FE.Estimate(1);
        subj.p1 = FE.pValue(2);
        subj.p0 = FE.pValue(1);

        if(subj.p1<0.05)
            a_sig_subjs{a_sig_cnt} = subj;
            a_sig_cnt = a_sig_cnt + 1;
        else
            a_non_subjs{a_non_cnt} = subj;
            a_non_cnt = a_non_cnt +1;
        end

        % ----------------------------------------
        cnt = cnt + 1;
        subjects = [subjects;repmat(cnt,length(gs_prds.v_cls_hd), ...
                                    1)];
    catch
        % do nothing
        logger(['  -Could not find/load prds for: ',subj_study,'_', ...
                name],proj.path.logfile);
    end

end

if(v_sig_cnt>1)
    save([proj.path.analysis.ex_gs_vs_gm,'v_sig_subjs.mat'],'v_sig_subjs');
end

if(v_non_cnt>1)
    save([proj.path.analysis.ex_gs_vs_gm,'v_sig_subjs.mat'],'v_non_subjs');
end

if(a_sig_cnt>1)
    save([proj.path.analysis.ex_gs_vs_gm,'v_sig_subjs.mat'],'a_sig_subjs');
end

if(a_non_cnt>1)
    save([proj.path.analysis.ex_gs_vs_gm,'v_sig_subjs.mat'],'a_non_subjs');
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

%% Save out model
save([proj.path.analysis.ex_gs_vs_gm,'v_mdl.mat'],'mdl');

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
        proj.param.plot.light_grey);
hold on;

%% Individual plots
if(v_sig_cnt>1)
for i=1:numel(v_non_subjs)
    plot(v_non_subjs{i}.stim,v_non_subjs{i}.stim*v_non_subjs{i}.b1+ ...
         v_non_subjs{i}.b0,'Color',proj.param.plot.light_grey,'LineWidth',1);
end
end

if(v_non_cnt>1)
for i=1:numel(v_sig_subjs)
    plot(v_sig_subjs{i}.stim,v_sig_subjs{i}.stim*v_sig_subjs{i}.b1+ ...
         v_sig_subjs{i}.b0,'Color',proj.param.plot.dark_grey,'LineWidth',2);
end
end

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

%% Save out model
save([proj.path.analysis.ex_gs_vs_gm,'a_mdl.mat'],'mdl');

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
        proj.param.plot.light_grey);
hold on;

%% Individual plots
if(a_non_cnt>1)
    for i=1:numel(a_non_subjs)
        plot(a_non_subjs{i}.stim,a_non_subjs{i}.stim*a_non_subjs{i}.b1+ ...
             a_non_subjs{i}.b0,'Color',proj.param.plot.light_grey,'LineWidth',1);
    end
end

if(a_sig_cnt>1)
    for i=1:numel(a_sig_subjs)
        plot(a_sig_subjs{i}.stim,a_sig_subjs{i}.stim*a_sig_subjs{i}.b1+ ...
             a_sig_subjs{i}.b0,'Color',proj.param.plot.dark_grey,'LineWidth',2);
    end
end

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

%% clean up
close all;

