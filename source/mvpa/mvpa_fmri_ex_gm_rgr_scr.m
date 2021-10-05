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
logger(['Intra-subject LOOCV MVPA RGR GM Features -> SCRs'],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.mvpa.fmri_ex_gm_rgr_scr]);
    eval(['! rm -rf ',proj.path.mvpa.fmri_ex_gm_rgr_scr]);
    disp(['Creating ',proj.path.mvpa.fmri_ex_gm_rgr_scr]);
    eval(['! mkdir ',proj.path.mvpa.fmri_ex_gm_rgr_scr]);
end

%% ----------------------------------------
%% Load labels;
label_id = load([proj.path.trg.ex,'stim_ids.txt']);

%% ----------------------------------------
%% load subjs
subjs = proj.process.subjs;

%% ----------------------------------------
%% iterate over study subjects

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
        all_img = base_img(in_brain,:)';
        
        %% Subselect extrinsic data
        ex_id = find(label_id==proj.param.trg.ex_id);
        ex_img = all_img(ex_id,:);
        
        %% Peform quality check of generated features
        qlty = check_gm_img_qlty(ex_img);
        
        if(qlty.ok)
            
            %% ----------------------------------------
            %% Train SCR
            
            %% Load SCR values
            load([proj.path.betas.scr_ex_beta,subj_study,'_',name,'_ex_betas.mat']);
            
            %%Change name to handle missing SCR
            scr_ex_img = ex_img;
            scr_ex_betas = [ex_betas.id1';ex_betas.id2'];
            
            %%Adjust training data to handle mising SCR
            if(isempty(ex_betas.id1))
                scr_ex_img = ex_img(46:90,:); 
            end
            
            if(isempty(ex_betas.id2))
                scr_ex_img = ex_img(1:45,:);  
            end

            %% Fit model
            [out,trg,mdl,stats] = regress_intra_loocv(zscore(scr_ex_img),zscore(scr_ex_betas),...
                                                      proj.param.mvpa.kernel);

            measures = [measures;trg];
            predictors = [predictors;out];
            subjects = [subjects;repmat(i,numel(out),1)];
            
            %% Fit model
            scr_model = fitrsvm(zscore(scr_ex_img),zscore(scr_ex_betas),...
                                'KernelFunction',proj.param.mvpa.kernel);
            
            %% Save model
            save([proj.path.mvpa.fmri_ex_gm_rgr_scr,subj_study,'_', ...
                  name,'_ex_gm_rgr_scr_model.mat'],'scr_model');

            % ----------------------------------------
            % Quality control below
            
            % build individual subject structures
            subj = struct();
            subj.study = subj_study;
            subj.name = name;
            
            [b stat] = robustfit(out,trg);
            subj.stim = out;
            subj.b1 = b(2); % slope
            subj.b0 = b(1); % intercept
            subj.p1 = stat.p(2); %slope
            subj.p0 = stat.p(1); %intercept

            % sort subjects by significance
            if(subj.p1<0.05)
                sig_subjs{sig_cnt} = subj;
                sig_cnt = sig_cnt + 1;
            else
                non_subjs{non_cnt} = subj;
                non_cnt = non_cnt + 1;
            end

        end
        
    catch
        disp(['   MVPA Error: possible missing beta series']);
    end

end


%% ----------------------------------------
%% save out subject groups

if(exist('sig_subjs'))
    save([proj.path.mvpa.fmri_ex_gm_rgr_scr,'sig_subjs.mat'],'sig_subjs');
end

if(exist('non_subjs'))
    save([proj.path.mvpa.fmri_ex_gm_rgr_scr,'non_subjs.mat'],'non_subjs');
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
    logger('   random effects matter',proj.path.logfile);
    mdl = mdl_re;
else
    logger('   random effects DO NOT matter',proj.path.logfile);
end
logger(' ',proj.path.logfile);

% save out model
save([proj.path.mvpa.fmri_ex_gm_rgr_scr,'rgr_scr_mdl.mat'],'mdl');

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
%% Format figure
ymin = -3;
ymax = 3;
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

xlabel('Predicted SCR Beta Coefficients');
ylabel('SCR Beta Coefficients');

%% ----------------------------------------
%% explot hi-resolution figure
export_fig 'EX_mvpa_predicted_scr_summary.png' -r300  
eval(['! mv ',proj.path.code,'EX_mvpa_predicted_scr_summary.png ',proj.path.fig]);

%% clean-up
close all;
