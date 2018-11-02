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
logger(['*************************************************'], ...
       proj.path.logfile);
logger(['Analyzing SCR responses to EX stimuli            '], ...
       proj.path.logfile);
logger(['*************************************************'], ...
       proj.path.logfile);

%% plot parameters
axisLabelFontSize = 18;
circleSize = 10;
white = [1,1,1];
light_grey = [.8,.8,.8];
dark_grey = [.6,.6,.6];

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% ----------------------------------------
%% Load labels;
v_label = load([proj.path.trg.ex,'stim_v_labs.txt']);
a_label = load([proj.path.trg.ex,'stim_a_labs.txt']);
label_id = load([proj.path.trg.ex,'stim_ids.txt']);
v_score = load([proj.path.trg.ex,'stim_v_scores.txt']);
a_score = load([proj.path.trg.ex,'stim_a_scores.txt']);

%% Adjust for extrinsic presentations
v_score = v_score(find(label_id==proj.param.trg.ex_id));
a_score = a_score(find(label_id==proj.param.trg.ex_id));

%% ----------------------------------------
%% scatter the underlying stim and feel
measures = [];
predictors = [];
subjects = [];

for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    % debug
    logger([subj_study,'_',name],proj.path.logfile);

    try
        load([proj.path.betas.scr_beta,subj_study,'_',name,'_ex_betas.mat']);
    catch
        logger('    Could not find scr beta file for processing.', ...
               proj.path.logfile);
    end

    scr_betas = [ex_betas.id1,ex_betas.id2];
    scr_a_score = a_score;

    %% ----------------------------------------
    %% Control for missing Identify runs
    if(isempty(ex_betas.id1))
        scr_a_score = a_score(46:90);
    end

    if(isempty(ex_betas.id2))
        scr_a_score = a_score(1:45);
    end

    %% ****************************************
    %% Remove hardcoding of the indices covered
    %% by runs 1 and 2 of the extrinsic stimuli
    %%
    %% TICKET
    %% ****************************************

    if(~isempty(scr_betas))

        measures = [measures;scr_a_score];
        predictors = [predictors;scr_betas'];
        subjects = [subjects;repmat(i,numel(scr_betas),1)];

    end

end



%% ----------------------------------------
%% Group GLMM fit
measures = double(measures);
predictors = double(predictors-mean(predictors));
subjects = double(subjects);

tbl = table(measures,predictors,subjects,'VariableNames',{'trg', ...
                    'pred','subj'});

mdl_fe = fitlme(tbl,['trg ~ 1 + pred']);
mdl_re= fitlme(tbl,['trg ~ 1 + pred + (1+pred|subj)']);

%%Explore random effects across model types
fe_v_re = compare(mdl_fe,mdl_re);

mdl = mdl_fe;
if(fe_v_re.pValue<0.05);
    disp('   random effects matter');
    mdl = mdl_re;
else
    disp('   random effects DO NOT matter');
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
scatter(predictors,measures,10,'MarkerFaceColor', ...
        proj.param.plot.white,'MarkerEdgeColor', ...
        proj.param.plot.very_light_grey);
hold on;

%% ----------------------------------------
%% format figure
ymin = 1;
ymax = 9;
xmin = -3;
xmax = 3;

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

xlabel('SCR Beta Coefficients');
ylabel('Normative Arousal Scores');

%% ----------------------------------------
%% explot hi-resolution figure
export_fig 'EX_scr_summary.png' -r300  
eval(['! mv ',proj.path.code,'EX_scr_summary.png ',proj.path.fig]);
