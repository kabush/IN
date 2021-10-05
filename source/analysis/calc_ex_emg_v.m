%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

function calc_ex_emg_v(proj,var_name,fstd)

%% ----------------------------------------
%% load subjs
subjs = proj.process.subjs; 

%% ----------------------------------------
%% Load labels;
label_id = load([proj.path.trg.ex,'stim_ids.txt']);
v_score = load([proj.path.trg.ex,'stim_v_scores.txt']);

%% Adjust for extrinsic presentations
v_score = v_score(find(label_id==proj.param.trg.ex_id));

%% Adjust for thresholding
v_score = v_score-mean(v_score);
v_std = std(v_score);

%% ----------------------------------------
%% scatter the underlying stim and feel
measures = [];
predictors = [];
subjects = [];

sig_cnt = 1;
non_cnt = 1;

for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    % debug
    logger([subj_study,'_',name],proj.path.logfile);

    try
        load([proj.path.betas.emg_ex_beta,subj_study,'_',name,'_ex_betas.mat']);
    catch
        logger('    Could not find emg beta file for processing.', ...
               proj.path.logfile);
    end

    emg_betas = [eval(['ex_betas.',var_name,'.id1']);...
                 eval(['ex_betas.',var_name,'.id2'])];

    
    keep_ids = find(abs(v_score)>fstd*v_std);
    emg_betas = zscore(emg_betas(keep_ids));
    emg_v_score = zscore(v_score(keep_ids));

    %% ****************************************
    %% Remove hardcoding of the indices covered
    %% by runs 1 and 2 of the extrinsic stimuli
    %% ****************************************
    if(~isempty(emg_betas))

        measures = [measures;emg_v_score];
        predictors = [predictors;zscore(emg_betas)];
        subjects = [subjects;repmat(i,numel(emg_betas),1)];

        % ----------------------------------------
        % Quality control below
        
        % build individual subject structures
        subj = struct();
        subj.study = subj_study;
        subj.name = name;
        
        [b stat] = robustfit(emg_betas,emg_v_score);
        subj.stim = zscore(emg_betas);
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
        
    end

end

%% ----------------------------------------
%% save out subject groups

if(exist('sig_subjs'))
    save([proj.path.analysis.ex_emg_v,'sig_subjs.mat'],'sig_subjs');
end

if(exist('non_subjs'))
    save([proj.path.analysis.ex_emg_v,'non_subjs.mat'],'non_subjs');
end

%% ----------------------------------------
%% Group GLMM fit
measures = double(measures);
predictors = double(predictors-mean(predictors));
subjects = double(subjects);

size(predictors)
size(measures)
size(subjects)

tbl = table(measures,predictors,subjects,'VariableNames',{'trg', ...
                    'pred','subj'});

mdl_fe = fitlme(tbl,['trg ~ 1 + pred']);
mdl_re= fitlme(tbl,['trg ~ 1 + pred + (1+pred|subj)']);

%%Explore random effects across model types
fe_v_re = compare(mdl_fe,mdl_re);

mdl = mdl_fe;
if(fe_v_re.pValue<0.05);
    logger('   random effects matter',proj.path.logfile);
    mdl = mdl_re;
else
    logger('   random effects DO NOT matter',proj.path.logfile);
end
logger(' ',proj.path.logfile);

%save out the model
save([proj.path.analysis.ex_emg_v,'ex_emg_v_',var_name,'_mdl.mat'],'mdl');

%% ----------------------------------------
%% Examine Main Effect
[~,~,FE] = fixedEffects(mdl);
if(FE.pValue(2)<0.05)
    disp('Fixed Effects are significant');
    disp(['  p=',num2str(FE.pValue(2))]);
end

%% ----------------------------------------
%% compute effect size
Rsqr = mdl.Rsquared.Adjusted;
Fsqr = Rsqr/(1-Rsqr);
logger(['  Rsqr=',num2str(Rsqr)],proj.path.logfile);
logger(['  Fsqr=',num2str(Fsqr)],proj.path.logfile);
disp(' ');

figure(1)
set(gcf,'color','w');

%% ----------------------------------------
%% plot all the datapoints
scatter(predictors,measures,10,'MarkerFaceColor', ...
        proj.param.plot.white,'MarkerEdgeColor', ...
        proj.param.plot.light_grey);
hold on;

%% ----------------------------------------
%% overlay the individual plots
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
%% format figure
ymin = -2;
ymax = 2;
xmin = -2.5;
xmax = 2.5;

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

% xlabel('EMG Beta Coefficients');
% ylabel('Normative Valence Scores');

%% ----------------------------------------
%% explot hi-resolution figure
export_fig 'EX_emg_summary.png' -r300  
eval(['! mv ',proj.path.code,'EX_emg_summary.png ',proj.path.fig,...
     'EX_emg_',var_name,'_summary.png']);

% clean up
close all;
