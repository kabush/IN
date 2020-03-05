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
    disp(['Removing ',proj.path.analysis.in_emg]);
    eval(['! rm -rf ',proj.path.analysis.in_emg]);
    disp(['Creating ',proj.path.analysis.in_emg]);
    eval(['! mkdir ',proj.path.analysis.in_emg]);
end

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
trajs = [];

for i = 1:numel(subjs) %ONLY CTM subjs have EMG

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    % debug
    logger([subj_study,'_',name],proj.path.logfile);

    try
        load([proj.path.betas.emg_in_beta,subj_study,'_',name,'_in_betas.mat']);
        load([proj.path.betas.emg_in_beta,subj_study,'_',name,'_feel_betas.mat']);

        % Construct EMG of stimulus
        in_zygo_betas = zscore([in_betas.zygo.id1;in_betas.zygo.id2]);
        in_corr_betas = zscore([in_betas.corr.id1;in_betas.corr.id2]);
        
        % Construct EMG of VR
        feel_zygo_betas = [feel_betas.zygo.id1;feel_betas.zygo.id2];
        feel_corr_betas = [feel_betas.corr.id1;feel_betas.corr.id2];
        zygo_measures = [zygo_measures;zscore(reshape(feel_zygo_betas',prod(size(feel_zygo_betas)),1))];
        corr_measures = [corr_measures;zscore(reshape(feel_corr_betas',prod(size(feel_corr_betas)),1))];

        % ----------------------------------------
        % Reshape for use in GLMM 

        % zygo predictors
        in_zygo_betas_box = repmat(in_zygo_betas,1,4);
        in_zygo_preds = reshape(in_zygo_betas_box',prod(size(in_zygo_betas_box)),1);
        zygo_predictors = [zygo_predictors;in_zygo_preds];
        
        % corr predictors
        in_corr_betas_box = repmat(in_corr_betas,1,4);
        in_corr_preds = reshape(in_corr_betas_box',prod(size(in_corr_betas_box)),1);
        corr_predictors = [corr_predictors;in_corr_preds];
        
        % controls
        subjects = [subjects;repmat(i,4*numel(in_zygo_betas),1)];
        traj_box = repmat(1:4,size(in_zygo_betas,1),1);
        traj_rshp = reshape(traj_box',prod(size(traj_box)),1);
        trajs = [trajs;traj_rshp];

    catch
        logger('    Could not find emg beta file for processing.',proj.path.logfile);
    end

end

logger('----------------------------------------',proj.path.logfile);
logger(' ZYGOMATICUS ANALYSIS                   ',proj.path.logfile);
logger('----------------------------------------',proj.path.logfile);

%% ----------------------------------------
%% Group GLMM fit
measures = double(zygo_measures);
predictors = double(zygo_predictors);
subjects = double(subjects);
trajs = double(trajs);

tbl = table(measures,predictors,trajs,subjects,'VariableNames',{'trg', ...
                    'pred','traj','subj'});

mdl_fe = fitlme(tbl,['trg ~ 1 + traj + pred']);
mdl_re= fitlme(tbl,['trg ~ 1 + pred + traj + (1+traj|subj) + (1+pred|subj)']);

%%Explore random effects across model types
fe_v_re = compare(mdl_fe,mdl_re);

mdl = mdl_fe;
logger('Testing Random Effects',proj.path.logfile);
if(fe_v_re.pValue<0.05);
    logger('  random effects matter',proj.path.logfile);
    mdl = mdl_re;
else
    logger('  random effects DO NOT matter',proj.path.logfile);
end
logger(' ',proj.path.logfile);

% save out the fit
save([proj.path.analysis.in_emg,'zygo_mdl.mat'],'mdl');

%% ----------------------------------------
%% Examine Main Effect
[~,~,FE] = fixedEffects(mdl);
if(FE.pValue(2)<0.05)
    logger('Fixed Effects are significant',proj.path.logfile);
    logger(['  p=',num2str(FE.pValue(2))],proj.path.logfile);
end

%% ----------------------------------------
%% compute effect size
Rsqr = mdl.Rsquared.Adjusted;
Fsqr = Rsqr/(1-Rsqr);
logger(['  Rsqr=',num2str(Rsqr)],proj.path.logfile);
logger(['  Fsqr=',num2str(Fsqr)],proj.path.logfile);
logger(' ',proj.path.logfile);

figure(1)
set(gcf,'color','w');

%% ----------------------------------------
%% plot all the datapoints
scatter(predictors,measures,10,'MarkerFaceColor', ...
        proj.param.plot.white,'MarkerEdgeColor', ...
        proj.param.plot.light_grey);
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
plot(vseq,y_hat,'r-','LineWidth',3);

xlim([xmin,xmax]);
ylim([ymin,ymax]);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

%% ----------------------------------------
%% explot hi-resolution figure
export_fig 'IN_emg_zygo_summary.png' -r300  
eval(['! mv ',proj.path.code,'IN_emg_zygo_summary.png ',proj.path.fig]);

%% clean up
close all;

logger('----------------------------------------',proj.path.logfile);
logger(' CORRUGATOR ANALYSIS                    ',proj.path.logfile);
logger('----------------------------------------',proj.path.logfile);

%% ----------------------------------------
%% Group GLMM fit
measures = double(corr_measures);
predictors = double(corr_predictors);
subjects = double(subjects);
trajs = double(trajs);

tbl = table(measures,predictors,trajs,subjects,'VariableNames',{'trg', ...
                    'pred','traj','subj'});

mdl_fe = fitlme(tbl,['trg ~ 1 + traj + pred']);
mdl_re= fitlme(tbl,['trg ~ 1 + pred + traj + (1+traj|subj) + (1+pred|subj)']);

%%Explore random effects across model types
fe_v_re = compare(mdl_fe,mdl_re);

mdl = mdl_fe;
logger('Testing Random Effects',proj.path.logfile);
if(fe_v_re.pValue<0.05);
    logger('  random effects matter',proj.path.logfile);
    mdl = mdl_re;
else
    logger('  random effects DO NOT matter',proj.path.logfile);
end
logger(' ',proj.path.logfile);

% save out the fit
save([proj.path.analysis.in_emg,'corr_mdl.mat'],'mdl');

%% ----------------------------------------
%% Examine Main Effect
[~,~,FE] = fixedEffects(mdl);
if(FE.pValue(2)<0.05)
    logger('Fixed Effects are significant',proj.path.logfile);
    logger(['  p=',num2str(FE.pValue(2))],proj.path.logfile);
end

%% ----------------------------------------
%% compute effect size
Rsqr = mdl.Rsquared.Adjusted;
Fsqr = Rsqr/(1-Rsqr);
logger(['  Rsqr=',num2str(Rsqr)],proj.path.logfile);
logger(['  Fsqr=',num2str(Fsqr)],proj.path.logfile);
logger(' ',proj.path.logfile);

figure(2)
set(gcf,'color','w');

%% ----------------------------------------
%% plot all the datapoints
scatter(predictors,measures,10,'MarkerFaceColor', ...
        proj.param.plot.white,'MarkerEdgeColor', ...
        proj.param.plot.light_grey);
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
plot(vseq,y_hat,'r-','LineWidth',3);

xlim([xmin,xmax]);
ylim([ymin,ymax]);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

%% ----------------------------------------
%% explot hi-resolution figure
export_fig 'IN_emg_corr_summary.png' -r300  
eval(['! mv ',proj.path.code,'IN_emg_corr_summary.png ', ...
      proj.path.fig]);

%% clean up
close all;


% %% ----------------------------------------
% %% ----------------------------------------
% %% Compare REST vs IN
% 
% rest_b = [];
% in_b = [];
% 
% % gather REST fits
% if(exist('non_subjs'))
%     for i =1:numel(non_subjs)
%         rest_b = [rest_b,non_subjs{i}.b1];
%     end
% end
% 
% if(exist('sig_subjs'))
%     for i =1:numel(sig_subjs)
%         rest_b = [rest_b,sig_subjs{i}.b1];
%     end
% end
% 
% % clean up
% clear non_subjs;
% clear sig_subjs;
% 
% % gather IN fits
% load([proj.path.analysis.vr_skill,var_name,'_sig_subjs.mat']);
% load([proj.path.analysis.vr_skill,var_name,'_non_subjs.mat']);
% 
% if(exist('non_subjs'))
%     for i =1:numel(non_subjs)
%         in_b = [in_b,non_subjs{i}.b1];
%     end
% end
% 
% if(exist('sig_subjs'))
%     for i =1:numel(sig_subjs)
%         in_b = [in_b,sig_subjs{i}.b1];
%     end
% end
% 
% % analysis
% logger(['*************************************'],proj.path.logfile);
% logger(['Comparison of IN ctrl to REST entrain'],proj.path.logfile);
% logger(['  Median IN beta: ',num2str(median(in_b))],proj.path.logfile);
% logger(['  Median RST beta: ',num2str(median(rest_b))],proj.path.logfile);
% [p h] = ranksum(in_b,rest_b);
% logger(['  Ranksum(IN_beta vs RST_beta), p=',num2str(p)],proj.path.logfile);
% 
% % plot effect
% 
% figure(1)
% set(gcf,'color','w');
% hold on;
% h1 = histogram(in_b,'BinWidth',0.1); 
% h2 = histogram(rest_b,'BinWidth',0.1,'FaceColor','red');
% hold off;
% 
% fig = gcf;
% ax = fig.CurrentAxes;
% ax.FontSize = proj.param.plot.axisLabelFontSize;
% 
% export_fig entrain_hist.png -r300  
% eval(['! mv ',proj.path.code,'entrain_hist.png ',proj.path.fig,'REST_vs_IN_',var_name,'_entrain_hist.png']);
% 
% % clean-up
% close all;



