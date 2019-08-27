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
logger(['*************************************************'],proj.path.logfile);
logger(['Computing IN Cog Control EVC Models                  '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%%%% %% ----------------------------------------
%%%% %% Set-up Directory Structure for fMRI betas
%%%% if(proj.flag.clean_build)
%%%%     disp(['Removing ',proj.path.ctrl.in_acc_activ]);
%%%%     eval(['! rm -rf ',proj.path.ctrl.in_acc_activ]);
%%%%     disp(['Creating ',proj.path.ctrl.in_acc_activ]);
%%%%     eval(['! mkdir ',proj.path.ctrl.in_acc_activ]);
%%%% end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% Storage for analysis
all_b_err = [];
all_b_cnf = [];
all_b_pel = [];
all_b_pro = [];
all_b_evc = [];

measures = [];
% predictors = [];
errs = [];
cnfs = [];
pels = [];
pros = [];
evcs = [];
subjects = [];

subj_cnt = 0;

%% ----------------------------------------
%% Transform beta-series into affect series {v,a}
for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    % log processing of subject
    logger([subj_study,'_',name],proj.path.logfile);

    data_exist = 0;
    try

        % load dynamics
        load([proj.path.ctrl.in_err_mdl,subj_study,'_',name,'_mdls.mat']);
        
        %% Data is present
        data_exist = 1;
        
    catch
        logger(['   -predictions do not exist'],proj.path.logfile);
    end

    if(data_exist)

        subj_cnt = subj_cnt+1;

        %% ----------------------------------------
        %% Find GM intersection with ACC ICA

        % Load GM mask
        gm_nii = load_nii([proj.path.mri.gm_mask,subj_study,'.',name,'.gm.nii']);
        mask = double(gm_nii.img);
        brain_size=size(mask);
        mask = reshape(mask,brain_size(1)*brain_size(2)*brain_size(3),1);
        in_brain=find(mask==1);  

        % Load ACCI ICA mask
        ica_nii = load_untouch_nii([proj.path.ctrl.in_ica,'sng_orient_thresh_zstatd70_17_3x3x3.nii.gz']);
        ica = double(ica_nii.img);
        brain_size=size(ica);
        ica = reshape(ica,brain_size(1)*brain_size(2)*brain_size(3),1);
        in_brain_this_ica = find(ica>0);

        % Intersection
        in_brain_this_ica = intersect(in_brain,in_brain_this_ica);
        disp(['  icaACC|GM voxels: ',num2str(numel(unique(in_brain_this_ica)))]);

        % Load beta-series
        base_nii = load_untouch_nii([proj.path.betas.fmri_in_beta,subj_study,'_',name,'_lss.nii']);
        brain_size = size(base_nii.img);
        
        % Vectorize the base image
        base_img = vec_img_2d_nii(base_nii);
        base_img = reshape(base_img,brain_size(1)*brain_size(2)*brain_size(3),brain_size(4))';

        %% ----------------------------------------
        %% Find ACC activation trajectory
        base_acc = mean(base_img(:,in_brain_this_ica),2);


        %% ----------------------------------------
        %% Load computational models

        %error
        load([proj.path.ctrl.in_err_mdl,subj_study,'_',name,'_mdls.mat']);
        err_mdls = mdls;

        %conflict
        load([proj.path.ctrl.in_cnf_mdl,subj_study,'_',name,'_mdls.mat']);
        cnf_mdls = mdls;

        %prediction err likelihood
        load([proj.path.ctrl.in_pel_mdl,subj_study,'_',name,'_mdls.mat']);
        pel_mdls = mdls;

        %prediction response outcome
        load([proj.path.ctrl.in_pro_mdl,subj_study,'_',name,'_mdls.mat']);
        pro_mdls = mdls;

        %expected value of control
        load([proj.path.ctrl.in_evc_mdl,subj_study,'_',name,'_mdls.mat']);
        evc_mdls = mdls;

        indx = reshape(err_mdls.a_indx',1, ...
                       prod(size(err_mdls.a_indx)));
        
        % true measure
        acc = base_acc(indx,1);

        % true predictors
        % err = reshape(abs(err_mdls.a_dcmp.err'),1,prod(size(err_mdls.a_dcmp.err)))'; %
        err = reshape(err_mdls.a_dcmp',1,prod(size(err_mdls.a_dcmp)))'; %***TICKET***
        cnf = reshape(cnf_mdls.a_dcmp',1,prod(size(cnf_mdls.a_dcmp)))';
        pel = reshape(pel_mdls.a_dcmp',1,prod(size(pel_mdls.a_dcmp)))';
        pro = reshape(pro_mdls.a_dcmp',1,prod(size(pro_mdls.a_dcmp)))';
        evc = reshape(evc_mdls.a_dcmp',1,prod(size(evc_mdls.a_dcmp)))';

        % true subject id
        subject = repmat(subj_cnt,numel(indx),1);

        % concatenate
        measures = [measures;zscore(acc)];
        errs = [errs;zscore(err)];
        cnfs = [cnfs;zscore(cnf)];
        pels = [pels;zscore(pel)];
        pros = [pros;zscore(pro)];
        evcs = [evcs;zscore(evc)];
        subjects = [subjects;subject];

    end
  
end

%% ----------------------------------------
%% Group GLMM fit
measures = double(zscore(measures));
errs = double(zscore(errs-mean(errs)));
cnfs = double(zscore(cnfs-mean(cnfs)));
pels = double(zscore(pels-mean(pels))); 
pros = double(zscore(pros-mean(pros))); 
evcs = double(zscore(evcs-mean(evcs))); 
subjects = double(subjects);

tbl = table(measures,errs,cnfs,pels,pros,evcs,subjects,'VariableNames',{'trg', ...
                    'err','cnf','pel','pro','evc','subj'});

mdl_fe = fitlme(tbl,['trg ~ 1 + err + cnf + pel + pro + evc']);
mdl_re= fitlme(tbl,['trg ~ 1 + err + cnf + pel + pro + evc + ' ...
                    '(err|subj) + (cnf|subj) + (pel|subj) + ' ...
                    '(pro|subj) + (evc|subj)']);

%%Explore random effects across model types
fe_a_re = compare(mdl_fe,mdl_re);

mdl = mdl_fe;
if(fe_a_re.pValue<0.05);
    disp('   random effects matter');
    mdl = mdl_re;
else
    disp('   random effects DO NOT matter');
end

disp(' ');

%% ----------------------------------------
%% Examine Main Effect
[~,~,FE] = fixedEffects(mdl);

for i =2:numel(FE.pValue)
    if(FE.pValue(i)<0.05)
        disp(['Fixed Effect [',num2str(i-1),'] is significant']);
        disp(['  p=',num2str(FE.pValue(i))]);
    end

end

%% ----------------------------------------
%% compute effect size
SS_res=sum((mdl.residuals).^2);
SS_tot=sum((measures-mean(measures)).^2);
Rsqr = 1-(SS_res/SS_tot);
Fsqr = Rsqr/(1-Rsqr);
logger(['Overall model fit'],proj.path.logfile);
logger(['  Rsqr=',num2str(Rsqr)],proj.path.logfile);
logger(['  Fsqr=',num2str(Fsqr)],proj.path.logfile);

disp(' ');


%% ----------------------------------------
%% format figure
ymin = -3;
ymax = 3;
xmin = -3;
xmax = 3; 
vseq = linspace(xmin,xmax);


%% ----------------------------------------
%% plot all the datapoints

%plot err
figure(1)
set(gcf,'color','w');
scatter(errs,measures,10,'MarkerFaceColor', ...
        proj.param.plot.white,'MarkerEdgeColor', ...
        proj.param.plot.very_light_grey);
hold on;
y_hat_err = FE.Estimate(1) + FE.Estimate(2)*vseq;
plot(vseq,y_hat_err,'r-','LineWidth',3);

xlim([xmin,xmax]);
ylim([ymin,ymax]);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

xlabel('Estimated Error');
ylabel('ACC activation');


%plot evc
figure(2)
set(gcf,'color','w');
scatter(evcs,measures,10,'MarkerFaceColor', ...
        proj.param.plot.white,'MarkerEdgeColor', ...
        proj.param.plot.very_light_grey);
hold on;
y_hat_evc = FE.Estimate(1) + FE.Estimate(6)*vseq;
plot(vseq,y_hat_evc,'r-','LineWidth',3);

xlim([xmin,xmax]);
ylim([ymin,ymax]);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

xlabel('Estimated Q-Value');
ylabel('ACC activation');



% %% ----------------------------------------
% %% explot hi-resolution figure
% export_fig 'R01_outcome_summary_a.png' -r300  
% eval(['! mv ',proj.path.code,'R01_outcome_summary_a.png ',proj.path.fig]);
