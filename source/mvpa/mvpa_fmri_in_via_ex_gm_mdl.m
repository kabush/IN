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
logger(['************************************************'],proj.path.logfile);
logger(['Validate EX modeling of IN data (Platt space) '],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.mvpa.fmri_in_via_ex_gm_mdl]);
    eval(['! rm -rf ',proj.path.mvpa.fmri_in_via_ex_gm_mdl]);
    disp(['Creating ',proj.path.mvpa.fmri_in_via_ex_gm_mdl]);
    eval(['! mkdir ',proj.path.mvpa.fmri_in_via_ex_gm_mdl]);
end

%% ----------------------------------------
%% Load labels;
v_label = load([proj.path.trg.ex,'stim_v_labs.txt']);
a_label = load([proj.path.trg.ex,'stim_a_labs.txt']);
label_id = load([proj.path.trg.ex,'stim_ids.txt']);
v_score = load([proj.path.trg.ex,'stim_v_scores.txt']);
a_score = load([proj.path.trg.ex,'stim_a_scores.txt']);

% % Scale scores to probabilities
v_score = v_score-1; 
v_score = v_score/9; 
a_score = a_score-1;
a_score = a_score/9; 

%% ----------------------------------------
%% load subjs
subjs = proj.process.subjs; 

all_v_trgs = [];
all_v_prds = [];
all_v_sbjs = [];

all_a_trgs = [];
all_a_prds = [];
all_a_sbjs = [];

clear v_sig_subjs;
clear v_non_subjs;

clear a_sig_subjs;
clear a_non_subjs;

v_sig_cnt = 1;
v_non_cnt = 1;
a_sig_cnt = 1;
a_non_cnt = 1;

%% ----------------------------------------
%% iterate over study subjects
for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    %% debug
    logger([subj_study,':',name],proj.path.logfile);

    try

        %% Load gray matter mask 
        gm_nii = load_nii([proj.path.mri.gm_mask,subj_study,'.',name,'.gm.nii']);
        mask = double(gm_nii.img);
        brain_size=size(mask);
        mask = reshape(mask,brain_size(1)*brain_size(2)*brain_size(3),1);
        in_brain=find(mask==1);  
        
        %% Load beta-series
        base_nii = load_nii([proj.path.betas.fmri_ex_beta,subj_study,'_',name,'_lss.nii']);
        brain_size = size(base_nii.img);
        
        %% Vectorize the base image
        base_img = vec_img_2d_nii(base_nii);
        base_img = reshape(base_img,brain_size(1)*brain_size(2)*brain_size(3),brain_size(4));
        
        %% Concatenate the MASKED base image
        subj_img = base_img(in_brain,:)';
        
        %% Concatenate all label/subj identifiers
        subj_id = repmat(id,numel(v_label),1);
        subj_i = repmat(i,numel(v_label),1);

        %% Subselect intrinsic  data
        in_id = find(label_id==proj.param.trg.in_id);

        in_id_good = in_id;
        in_img = subj_img(in_id_good,:);
        in_subj_id = subj_id(in_id_good,1);
        in_v_label = v_label(in_id_good,1);
        in_v_score = v_score(in_id_good,1);
        in_a_label = a_label(in_id_good,1);
        in_a_score = a_score(in_id_good,1);

        %% Peform quality check of generated features
        qlty = check_gm_img_qlty(in_img);
        
        if(qlty.ok)
            
            %% ----------------------------------------
            %% Predict IN VALENCE with EX model
            
            %% Load valence model
            load([proj.path.mvpa.fmri_ex_gm_mdl,subj_study,'_',name,'_v_model.mat']);

            %% Predict labels of IN stims
            [tst_predict,tst_score] = predict(v_model,in_img);
            prds = 1./(1+exp(-tst_score(:,2))); % Platt scale predictions

            %% Gather data and scale: within subject
            all_v_trgs = [all_v_trgs;in_v_score];
            all_v_prds = [all_v_prds;prds]; %zscore(tst_score(:,2))];
            all_v_sbjs = [all_v_sbjs;repmat(i,numel(in_v_score),1)];

            % ----------------------------------------
            % VALENCE Quality control below
            subj = struct();
            subj.study = subj_study;
            subj.name = name;
            
            measures = double(in_v_score);
            predictors = double(prds);
            
            tbl = table(measures,predictors,'VariableNames', ...
                        {'trgs','preds'});
             
            sbj_mdl_fe = fitlme(tbl,['trgs ~ 1 + preds']);
            sbj_mdl_fe.Rsquared.Adjusted
           
            % Extract Fixed effects
            [~,~,FE] = fixedEffects(sbj_mdl_fe);
            
            subj.stim = predictors;
            subj.b1 = FE.Estimate(2); % slope
            subj.b0 = FE.Estimate(1); % intercept
            subj.p1 = FE.pValue(2); %slope
            subj.p0 = FE.pValue(1); %intercept
            
            %Sort subjects by significance
            if(subj.p1<0.05)
                v_sig_subjs{v_sig_cnt} = subj;
                v_sig_cnt = v_sig_cnt + 1;
            else
                v_non_subjs{v_non_cnt} = subj;
                v_non_cnt = v_non_cnt + 1;
            end

            %% ----------------------------------------
            %% Predict IN AROUSAL with EX model
            
            %% Load valence model
            load([proj.path.mvpa.fmri_ex_gm_mdl,subj_study,'_',name,'_a_model.mat']);

            %% Predict labels of IN stims
            [tst_predict,tst_score] = predict(a_model,in_img);
            prds = 1./(1+exp(-tst_score(:,2))); % Platt scale predictions

            %% Gather data and scale: within subject
            all_a_trgs = [all_a_trgs;in_a_score];
            all_a_prds = [all_a_prds;prds]; %zscore(tst_score(:,2))];
            all_a_sbjs = [all_a_sbjs;repmat(i,numel(in_a_score),1)];

            % ----------------------------------------
            % AROUSAL Quality control below
            subj = struct();
            subj.study = subj_study;
            subj.name = name;
            
            measures = double(in_a_score);
            predictors = double(prds);
            
            tbl = table(measures,predictors,'VariableNames', ...
                        {'trgs','preds'});
             
            sbj_mdl_fe = fitlme(tbl,['trgs ~ 1 + preds']);
            sbj_mdl_fe.Rsquared.Adjusted
           
            % Extract Fixed effects
            [~,~,FE] = fixedEffects(sbj_mdl_fe);
            
            subj.stim = predictors;
            subj.b1 = FE.Estimate(2); % slope
            subj.b0 = FE.Estimate(1); % intercept
            subj.p1 = FE.pValue(2); %slope
            subj.p0 = FE.pValue(1); %intercept
            
            %Sort subjects by significance
            if(subj.p1<0.05)
                a_sig_subjs{a_sig_cnt} = subj;
                a_sig_cnt = a_sig_cnt + 1;
            else
                a_non_subjs{a_non_cnt} = subj;
                a_non_cnt = a_non_cnt + 1;
            end

        end
        
    catch
        logger(['  -MVPA Error: possible missing beta series'],proj.path.logfile);
    end

end

%% Save out sig and non-sig subjs
if(exist('v_sig_subjs'))
    save([proj.path.mvpa.fmri_in_via_ex_gm_mdl,'v_sig_subjs.mat'],'v_sig_subjs');
end
if(exist('v_non_subjs'))
    save([proj.path.mvpa.fmri_in_via_ex_gm_mdl,'v_non_subjs.mat'],'v_non_subjs');
end
if(exist('a_sig_subjs'))
    save([proj.path.mvpa.fmri_in_via_ex_gm_mdl,'a_sig_subjs.mat'],'a_sig_subjs');
end
if(exist('a_non_subjs'))
    save([proj.path.mvpa.fmri_in_via_ex_gm_mdl,'a_non_subjs.mat'],'a_non_subjs');
end


%% ----------------------------------------
%% Analyze VALENCE Effects

%% Group GLMM fit
measures = double(all_v_trgs);% -mean(all_v_trgs));
predictors = double(all_v_prds); %-mean(all_v_prds));
subjects = double(all_v_sbjs);

tbl = table(measures,predictors,subjects,'VariableNames',{'trg', ...
                    'pred','subj'});

mdl_fe = fitlme(tbl,['trg ~ 1 + pred']);
mdl_re= fitlme(tbl,['trg ~ 1 + pred + (1+pred|subj)']);

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

%% Examine Main Effect
[~,~,FE] = fixedEffects(mdl);
if(FE.pValue(2)<0.05)
    logger('Fixed Effects are significant',proj.path.logfile);
    logger(['  p=',num2str(FE.pValue(2))],proj.path.logfile);
end

%% compute effect size
Rsqr = mdl.Rsquared.Adjusted;
Fsqr = Rsqr/(1-Rsqr);
logger(['  Rsqr=',num2str(Rsqr)],proj.path.logfile);
logger(['  Fsqr=',num2str(Fsqr)],proj.path.logfile);
logger(' ',proj.path.logfile);

%% save out the fit
save([proj.path.mvpa.fmri_in_via_ex_gm_mdl,'v_model.mat'],'mdl');

figure(1)
set(gcf,'color','w');
hold on;

%% ----------------------------------------
%% plot all the datapoints
scatter(predictors,measures,10,'MarkerFaceColor', ...
        proj.param.plot.white,'MarkerEdgeColor', ...
        proj.param.plot.light_grey);

%% ----------------------------------------
%% overlay individual plots
if(exist('v_non_subjs'))
    for i =1:numel(v_non_subjs)
        plot(v_non_subjs{i}.stim,v_non_subjs{i}.stim*v_non_subjs{i}.b1+ ...
             v_non_subjs{i}.b0,'Color',proj.param.plot.light_grey,'LineWidth',1);
    end
end

if(exist('v_sig_subjs'))
    for i =1:numel(v_sig_subjs)
        plot(v_sig_subjs{i}.stim,v_sig_subjs{i}.stim*v_sig_subjs{i}.b1+ ...
             v_sig_subjs{i}.b0,'Color',proj.param.plot.dark_grey,'LineWidth',2);
    end
end

%% ----------------------------------------
%% format figure (probabilities)
ymin = 0;
ymax = .85;
xmin = 0;
xmax = 1;

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
export_fig 'IN_val_via_ex_mdl.png' -r300  
eval(['! mv ',proj.path.code,'IN_val_via_ex_mdl.png ',proj.path.fig]);


%% ----------------------------------------
%% Analyze AROUSAL Effects

%% Group GLMM fit
measures = double(all_a_trgs); %-mean(all_a_trgs));
predictors = double(all_a_prds); %-mean(all_a_prds));
subjects = double(all_a_sbjs);

tbl = table(measures,predictors,subjects,'VariableNames',{'trg', ...
                    'pred','subj'});

mdl_fe = fitlme(tbl,['trg ~ 1 + pred']);
mdl_re= fitlme(tbl,['trg ~ 1 + pred + (1+pred|subj)']);

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

%% Examine Main Effect
[~,~,FE] = fixedEffects(mdl);
if(FE.pValue(2)<0.05)
    logger('Fixed Effects are significant',proj.path.logfile);
    logger(['  p=',num2str(FE.pValue(2))],proj.path.logfile);
end

%% compute effect size
Rsqr = mdl.Rsquared.Adjusted;
Fsqr = Rsqr/(1-Rsqr);
logger(['  Rsqr=',num2str(Rsqr)],proj.path.logfile);
logger(['  Fsqr=',num2str(Fsqr)],proj.path.logfile);
logger(' ',proj.path.logfile);

%% save out the fit
save([proj.path.mvpa.fmri_in_via_ex_gm_mdl,'a_model.mat'],'mdl');

figure(2)
set(gcf,'color','w');
hold on;

%% ----------------------------------------
%% plot all the datapoints
scatter(predictors,measures,10,'MarkerFaceColor', ...
        proj.param.plot.white,'MarkerEdgeColor', ...
        proj.param.plot.light_grey);

%% ----------------------------------------
%% overlay individual plots
if(exist('a_non_subjs'))
    for i =1:numel(a_non_subjs)
        plot(a_non_subjs{i}.stim,a_non_subjs{i}.stim*a_non_subjs{i}.b1+ ...
             a_non_subjs{i}.b0,'Color',proj.param.plot.light_grey,'LineWidth',1);
    end
end

if(exist('a_sig_subjs'))
    for i =1:numel(a_sig_subjs)
        plot(a_sig_subjs{i}.stim,a_sig_subjs{i}.stim*a_sig_subjs{i}.b1+ ...
             a_sig_subjs{i}.b0,'Color',proj.param.plot.dark_grey,'LineWidth',2);
    end
end

%% ----------------------------------------
%% format figure (probabilities)
ymin = 0;
ymax = .85;
xmin = 0;
xmax = 1;

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
export_fig 'IN_aro_via_ex_mdl.png' -r300  
eval(['! mv ',proj.path.code,'IN_aro_via_ex_mdl.png ',proj.path.fig]);

%% clean up
close all;
