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

figure(1)
set(gcf,'color','w');

%% ----------------------------------------
%% scatter the underlying stim and feel
indv_b = [];

for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    % debug
    logger([subj_study,'_',name],proj.path.logfile);

    try
        load([proj.path.betas.scr_beta,subj_study,'_',name, ...
              '_ex_betas.mat']);
    catch
        logger('    Could not find scr beta file for processing.', ...
               proj.path.logfile);
    end

    scr_betas = [ex_betas.id1,ex_betas.id2];
    scr_a_score = a_score;

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

        %% scatter plot specific points        
        scatter(scr_betas,scr_a_score,10,'MarkerFaceColor',white, ...
                'MarkerEdgeColor',light_grey);
        hold on;
        
        %% robust fit
        [b stat] = robustfit(scr_betas,scr_a_score);
        indv_b = [indv_b;b'];
        disp(stat.p(2))
        
        plot(sort(scr_betas),sort(scr_betas)*b(2)+b(1),'r-');
        
    end

end

[h p ci stat] = ttest(indv_b(:,2));
disp(['***EX response beta ci=[',num2str(ci(1)),' ',num2str(ci(2)),['], ' ...
                    '' ...
                    'p='],num2str(p),'***']);
