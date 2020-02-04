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
    disp(['Removing ',proj.path.analysis.in_scr]);
    eval(['! rm -rf ',proj.path.analysis.in_scr]);
    disp(['Creating ',proj.path.analysis.in_scr]);
    eval(['! mkdir ',proj.path.analysis.in_scr]);
end

%% Initialize log section
logger(['*********************************************'],proj.path.logfile);
logger(['Analyzing SCR responses to IN stimuli        '],proj.path.logfile);
logger(['*********************************************'],proj.path.logfile);

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% ----------------------------------------
%% scatter the underlying stim and feel
measures = [];
predictors = [];
subjects = [];
trajs = [];

for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    % debug
    logger([subj_study,'_',name],proj.path.logfile);

    try
        load([proj.path.betas.scr_in_beta,subj_study,'_',name,'_in_betas.mat']);
        load([proj.path.betas.scr_in_beta,subj_study,'_',name,'_feel_betas.mat']);
    catch
        logger('    Could not find scr beta file for processing.',proj.path.logfile);
    end

    scr_in_betas = [in_betas.id1,in_betas.id2];
    scr_feel_betas = [feel_betas.id1;feel_betas.id2];
    
    % vectorize stimulus measures
    in_box = repmat(scr_in_betas',1,4);
    in_betas_tmp = reshape(in_box',prod(size(in_box)),1);

    % vectorize feel measures
    feel_betas_tmp = reshape(scr_feel_betas',prod(size(scr_feel_betas)),1);

    % build trajectory measure of time
    traj_box = repmat(1:4,numel(scr_in_betas),1);
    traj = reshape(traj_box',prod(size(traj_box)),1);


    if(~isempty(in_betas_tmp))

        measures = [measures;zscore(feel_betas_tmp)];
        predictors = [predictors;zscore(in_betas_tmp)];
        subjects = [subjects;repmat(i,numel(feel_betas_tmp),1)];
        trajs = [trajs;zscore(traj)];

    end

end

%% ----------------------------------------
%% Group GLMM fit
measures = double(measures);
predictors = double(predictors-mean(predictors));
trajs = double(trajs);
subjects = double(subjects);

tbl = table(measures,predictors,trajs,subjects,'VariableNames',{'trg', ...
                    'pred','traj','subj'});

mdl_fe = fitlme(tbl,['trg ~ 1 + pred + traj']);
mdl_re= fitlme(tbl,['trg ~ 1 + pred + traj + (1+pred|subj) + (1+traj|subj)']);

%%Explore random effects across model types
fe_v_re = compare(mdl_fe,mdl_re);

mdl = mdl_fe;
if(fe_v_re.pValue<0.05);
    logger('  random effects matter',proj.path.logfile);
    mdl = mdl_re;
else
    logger('  random effects DO NOT matter',proj.path.logfile);
end
logger(' ',proj.path.logfile);

% save out the fit
save([proj.path.analysis.in_scr,'scr_mdl.mat'],'mdl');

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
ymin = -3; %-2.88;
ymax = 3; %3.63;
xmin = -3; %-3.24;
xmax = 3; %3.61;

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

% xlabel('VR(cue) SCR responses');
% ylabel('VR(modulate) SCR responses');

%% ----------------------------------------
%% explot hi-resolution figure
export_fig 'IN_scr_summary.png' -r300  
eval(['! mv ',proj.path.code,'IN_scr_summary.png ',proj.path.fig]);

%% clean up
close all;