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

clear sig_subjs;
clear non_subjs;
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

        load([proj.path.betas.scr_in_beta,subj_study,'_',name,'_in_betas.mat']);
        load([proj.path.betas.scr_in_beta,subj_study,'_',name,'_feel_betas.mat']);
        
        % Construct EMG of stimulus
        in_scr = eval(['[in_betas.id1,in_betas.id2]']);;
        in_box = repmat(in_scr',1,4);
        feel_box = eval(['[feel_betas.id1;feel_betas.id2]']);
        traj_box = repmat(1:4,numel(in_scr),1);
        
        % Create symmetric boxes
        in_form = zscore(reshape(in_box',prod(size(in_box)),1));
        feel_form = zscore(reshape(feel_box',prod(size(feel_box)),1));
        traj_form = reshape(traj_box',prod(size(traj_box)),1);
        
        % Build data for group GLMM
        predictors = [predictors;in_form];
        measures = [measures;feel_form];
        subjects = [subjects;repmat(i,numel(in_form),1)];
        trajs = [trajs;traj_form];
        
        % ----------------------------------------
        % Quality control below
        
        % build individual subject structures
        subj = struct();
        subj.study = subj_study;
        subj.name = name;
        
        tbl = table(feel_form,in_form,traj_form,'VariableNames', ...
                    {'feel','stim','traj'});
        
        sbj_mdl_fe = fitlme(tbl,['feel ~ 1 + stim + traj']);
        
        %% Extract Fixed effects
        [~,~,FE] = fixedEffects(sbj_mdl_fe);
        
        subj.stim = in_form;
        subj.b1 = FE.Estimate(2); % slope
        subj.b0 = FE.Estimate(1); % intercept
        subj.p1 = FE.pValue(2); %slope
        subj.p0 = FE.pValue(1); %intercept
         
        %% sort subjects by significance
        if(subj.p1<0.05)
            sig_subjs{sig_cnt} = subj;
            sig_cnt = sig_cnt + 1
        else
            non_subjs{non_cnt} = subj;
            non_cnt = non_cnt + 1
        end
        
    catch
        logger('    Could not find scr beta file for processing.',proj.path.logfile);
    end
    
end

%% ----------------------------------------
%% save out subject groups
if(exist('sig_subjs'))
    save([proj.path.analysis.in_scr,'scr_sig_subjs.mat'],'sig_subjs');
end

if(exist('non_subjs'))
    save([proj.path.analysis.in_scr,'scr_non_subjs.mat'],'non_subjs');
end

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
save([proj.path.analysis.in_scr,'in_scr_mdl.mat'],'mdl');

%% ----------------------------------------
%% Plotting

figure(1)
set(gcf,'color','w');
hold on;

%% scatter raw data
scatter(predictors,measures,10,'MarkerFaceColor', ...
        proj.param.plot.white,'MarkerEdgeColor', ...
        proj.param.plot.light_grey);

%% ----------------------------------------
%% overlay the individual entrain plots
if(exist('non_subjs'))
    for i =1:numel(non_subjs)
        plot(non_subjs{i}.stim,non_subjs{i}.stim*non_subjs{i}.b1+ ...
             non_subjs{i}.b0,'Color',proj.param.plot.light_grey, ...
             'LineWidth',1);
    end
end

if(exist('sig_subjs'))
    for i =1:numel(sig_subjs)
        plot(sig_subjs{i}.stim,sig_subjs{i}.stim*sig_subjs{i}.b1+ ...
             sig_subjs{i}.b0,'Color',proj.param.plot.dark_grey, ...
             'LineWidth',2);
    end
end

%% ----------------------------------------
%% identify max/min x-range|y-rang
xmin = -3;
xmax = 3;
ymin = -3;
ymax = 3;

%% ----------------------------------------
%% overlay the group entrain plot
xseq = linspace(xmin,xmax);
y_hat = FE.Estimate(1) + FE.Estimate(2)*xseq;
plot(xseq,y_hat,'r-','LineWidth',3);

%% ----------------------------------------
%% format figure
xlim([xmin,xmax]);
ylim([ymin,ymax]);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

%% ----------------------------------------
%% export hi-resolution figure
export_fig in_scr_summary.png -r300  
eval(['! mv ',proj.path.code,'in_scr_summary.png ',proj.path.fig,'IN_scr_summary.png']);

% clean-up
close all;