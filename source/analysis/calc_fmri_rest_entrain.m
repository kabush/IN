%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

function calc_fmri_rest_entrain(proj,var_name)

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

figure(1)
set(gcf,'color','w');
hold on;

%% ----------------------------------------
%% scatter the underlying stim and feel
predictors = [];
trajs = [];
measures = [];
subjects = [];

sig_cnt = 1;
non_cnt = 1;

b_all = []; %% for power

for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    % log analysis of subject
    logger([subj_study,'_',name],proj.path.logfile);

    try

        path = [proj.path.mvpa.fmri_rest_via_ex_gm_mdl,subj_study,'_',name,'_slf_mu_',var_name,'_rest.mat'];

        %% Load REST trajectory
        load(path);
        data = eval(['slf_mu_',var_name,'_rest']);

        rest_start = proj.param.rest.n_trs_trans+1;
        rest_end = numel(data)-proj.param.rest.n_trs_tail-5; %shift
                                                             %left
                                                             %to
                                                             %allow
                                                             %for feel
        rest_ids = rest_start:rest_end;

        % plot(data(rest_start:rest_end));

        stim_ids = randsample(1:numel(rest_ids),proj.param.rest.n_resample)';
        stim_val = data(stim_ids);
        
        if(sum(isnan(data))==0)

            stim_box = repmat(stim_val,1,4);
            feel_box = [data(stim_ids+2),data(stim_ids+3),data(stim_ids+4),data(stim_ids+5)];
            traj_box = repmat(1:4,proj.param.rest.n_resample,1);
            
            %% reformat
            stim_form = reshape(stim_box',1,prod(size(stim_box)))';
            feel_form = reshape(feel_box',1,prod(size(feel_box)))';
            traj_form = reshape(traj_box',1,prod(size(traj_box)))';
            stim_form = zscore(stim_form);
            feel_form = zscore(feel_form);
            
            %% build data for group GLMM
            predictors = [predictors;stim_form]; 
            measures = [measures;feel_form];
            subjects = [subjects;repmat(i,numel(stim_form),1)];
            trajs = [trajs;traj_form];
            
            
            %% ----------------------------------------
            %% Quality control below
            
            %% build individual subject structures
            subj = struct();
            subj.study = subj_study;
            subj.name = name;
            
            tbl = table(feel_form,stim_form,traj_form,'VariableNames', ...
                        {'feel','stim','traj'});
            
            sbj_mdl_fe = fitlme(tbl,['feel ~ 1 + stim + traj']);
            
            %% Extract Fixed effects
            [~,~,FE] = fixedEffects(sbj_mdl_fe);
            
            subj.stim = zscore(stim_val)
            subj.b1 = FE.Estimate(2); % slope
            subj.b0 = FE.Estimate(1); % intercept
            subj.p1 = FE.pValue(2); %slope
            subj.p0 = FE.pValue(1); %intercept
            
            %% sort subjects by significance
            if(subj.p1<0.05)
                sig_subjs{sig_cnt} = subj;
                sig_cnt = sig_cnt + 1;
            else
                non_subjs{non_cnt} = subj;
                non_cnt = non_cnt + 1;
            end
        
        else
            logger(['  -Rst contains NaNs: ',subj_study,'_',name],proj.path.logfile);
        end
        
    catch
        % do nothing
        logger(['  -Could not find/load prds for: ',subj_study,'_', ...
                name],proj.path.logfile);
    end
        
end

%% ----------------------------------------
%% save out subject groups
if(exist('sig_subjs'))
    save([proj.path.analysis.fmri_rest_entrain,var_name,'_sig_subjs.mat'],'sig_subjs');
end

if(exist('non_subjs'))
    save([proj.path.analysis.fmri_rest_entrain,var_name,'_non_subjs.mat'],'non_subjs');
end


%% ----------------------------------------
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
save([proj.path.analysis.fmri_rest_entrain,var_name,'_rest_entrain_mdl.mat'],'mdl');

%% ----------------------------------------
%% ----------------------------------------
%% Plotting

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
xmin = -2.5;
xmax = 2.5;
ymin = -2;
ymax = 2;

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
export_fig entrain_summary.png -r300  
eval(['! mv ',proj.path.code,'entrain_summary.png ',proj.path.fig,'REST_',var_name,'_entrain_summary.png']);

% clean-up
close all;

%% ----------------------------------------
%% ----------------------------------------
%% Compare REST vs IN

rest_cnt = 1;

% gather REST fits
if(exist('non_subjs'))
    for i =1:numel(non_subjs)
        subj = struct();
        subj.study = non_subjs{i}.study;
        subj.name = non_subjs{i}.name;
        subj.b1 = non_subjs{i}.b1;
        rest_subjs{rest_cnt} = subj;
        rest_cnt = rest_cnt + 1;
    end
end

if(exist('sig_subjs'))
    for i =1:numel(sig_subjs)
        subj = struct();
        subj.study = sig_subjs{i}.study;
        subj.name = sig_subjs{i}.name;
        subj.b1 = sig_subjs{i}.b1;
        rest_subjs{rest_cnt} = subj;
        rest_cnt = rest_cnt + 1;
    end
end


% clean up
clear non_subjs;
clear sig_subjs;

% gather IN fits
load([proj.path.analysis.vr_skill,var_name,'_sig_subjs.mat']);
load([proj.path.analysis.vr_skill,var_name,'_non_subjs.mat']);

cmp_b = [];

if(exist('non_subjs'))
    for i =1:numel(non_subjs)
        for j=1:numel(rest_subjs)
            % disp([non_subjs{i}.study,'_',non_subjs{i}.name,':',rest_subjs{j}.study,'_',rest_subjs{j}.name]);
            if(strcmp(rest_subjs{j}.study,non_subjs{i}.study) & strcmp(rest_subjs{j}.name,non_subjs{i}.name))
                % disp(['   match**********']);
                cmp_b = [cmp_b,rest_subjs{j}.b1-non_subjs{i}.b1];
            end
        end
    end
end

if(exist('sig_subjs'))
    for i =1:numel(sig_subjs)
        for j=1:numel(rest_subjs)
            %     disp([sig_subjs{i}.study,'_',sig_subjs{i}.name,':',rest_subjs{j}.study,'_',rest_subjs{j}.name]);
            if(strcmp(rest_subjs{j}.study,sig_subjs{i}.study) & strcmp(rest_subjs{j}.name,sig_subjs{i}.name))
                %   disp(['   match**********']);
                cmp_b = [cmp_b,sig_subjs{i}.b1-rest_subjs{j}.b1];
            end
        end
    end
end

% analysis
logger(['*************************************'],proj.path.logfile);
logger(['Comparison of IN ctrl to REST entrain'],proj.path.logfile);
logger(['  Median cmp_beta: ',num2str(median(cmp_b))],proj.path.logfile);
[p,h,stats] = signrank(cmp_b);
logger(['  Signrank cmp_beta, p=',num2str(p)],proj.path.logfile);
r = stats.zval/sqrt(numel(cmp_b));
logger(['  ***EFFECT*** r(z/sqrt(n))=',num2str(r)],proj.path.logfile);


% plot effect
figure(1)
set(gcf,'color','w');
hold on;

colors = [1.0,0.0,0.0;
          0.8,0.8,0.8];

x = [cmp_b'];
g = repmat(1,numel(cmp_b),1);
g_labels = {' '};
kabBoxPlot(x,g,g_labels,colors,0.4)

fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = 22; %proj.param.plot.axisLabelFontSize;
set(gcf,'position',[0,0,300,600]);

export_fig entrain_hist.png -r300  
eval(['! mv ',proj.path.code,'entrain_hist.png ',proj.path.fig,'REST_vs_IN_',var_name,'_entrain_hist.png']);

% clean-up
close all;
