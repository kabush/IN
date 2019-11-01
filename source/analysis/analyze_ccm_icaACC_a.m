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
measures = [];
errs = [];
cnfs = [];
pels = [];
pros = [];
evcs = [];
trajs = [];
subjects = [];

%% ----------------------------------------
%% Transform beta-series into affect series {v,a}
subj_cnt = 0;
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
        ica_nii = load_untouch_nii([proj.path.ctrl.in_ica,'clst_sng_orient_thresh_zstatd70_17_3x3x3.nii.gz']);
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
        err = reshape(sqrt((err_mdls.a_dcmp').^2),1,prod(size(err_mdls.a_dcmp)))'; 
        cnf = reshape(cnf_mdls.a_dcmp',1,prod(size(cnf_mdls.a_dcmp)))';
        pel = reshape(pel_mdls.a_dcmp',1,prod(size(pel_mdls.a_dcmp)))';
        pro = reshape(pro_mdls.a_dcmp',1,prod(size(pro_mdls.a_dcmp)))';
        evc = reshape(evc_mdls.a_dcmp',1,prod(size(evc_mdls.a_dcmp)))';

        traj_box = repmat([1:4],30,1);
        traj = reshape(traj_box',1,prod(size(traj_box)))';

        % true subject id
        subject = repmat(subj_cnt,numel(indx),1);

        % concatenate
        measures = [measures;zscore(acc)];
        errs = [errs;zscore(err)];
        cnfs = [cnfs;zscore(cnf)];
        pels = [pels;zscore(pel)];
        pros = [pros;zscore(pro)];
        evcs = [evcs;zscore(evc)];
        trajs = [trajs;traj];
        subjects = [subjects;subject];

    end
  
end

%% ----------------------------------------
%% Assemble Measures
% measures = double(zscore(measures));
% errs = double(zscore(errs-mean(errs)));
% cnfs = double(zscore(cnfs-mean(cnfs)));
% pels = double(zscore(pels-mean(pels))); 
% pros = double(zscore(pros-mean(pros))); 
% evcs = double(zscore(evcs-mean(evcs))); 
% trajs = double(trajs);
% subjects = double(subjects);

measures = double(zscore(measures));
errs = double((errs-mean(errs)));
cnfs = double((cnfs-mean(cnfs)));
pels = double((pels-mean(pels))); 
pros = double((pros-mean(pros))); 
evcs = double((evcs-mean(evcs))); 
trajs = double(trajs);
subjects = double(subjects);



%% Group GLMM fit
% % tbl = table(measures,errs,cnfs,pels,pros,evcs,subjects,'VariableNames',{'trg', ...
% %                     'err','cnf','pel','pro','evc','subj'});
% tbl = table(measures,errs,cnfs,pels,pros,evcs,trajs,subjects,'VariableNames',{'trg', ...
%                     'err','cnf','pel','pro','evc','traj','subj'});
% 
% 
% mdl_fe = fitlme(tbl,['trg ~ 1 + err + cnf + pel + pro + evc']);
% % mdl_re= fitlme(tbl,['trg ~ 1 + err + cnf + pel + pro + evc + ' ...
% %                     '(err|subj) + (cnf|subj) + (pel|subj) + ' ...
% %                     '(pro|subj) + (evc|subj)']);
% mdl_re= fitlme(tbl,['trg ~ 1 + err + cnf + pel + pro + evc + ' ...
%                   '(err|traj) + (cnf|traj) + (pel|traj) + ' ...
%                     '(pro|traj) + (evc|traj) + ' ...
%                   '(err|subj) + (cnf|subj) + (pel|subj) + ' ...
%                     '(pro|subj) + (evc|subj)']);

tbl = table(measures,trajs,errs,cnfs,pels,pros,evcs,subjects,'VariableNames',{'trg', ...
                    'traj','err','cnf','pel','pro','evc','subj'});

mdl_fe = fitlme(tbl,['trg ~ 1 + traj + err + cnf + pel + pro + evc']);
mdl_re= fitlme(tbl,['trg ~ 1 + traj + err + cnf + pel + pro + evc + ' ...
                    '(err|subj) + (cnf|subj) + (pel|subj) + ' ...
                    '(pro|subj) + (evc|subj)']);





disp(' ');

%%Explore random effects across model types
fe_a_re = compare(mdl_fe,mdl_re);

mdl = mdl_fe;
if(fe_a_re.pValue<0.05);
    disp('Random effects matter');
    mdl = mdl_re;
else
    disp('Random effects DO NOT matter');
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
Rsqr = mdl.Rsquared.Ordinary; ;
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


%%%%% DEBUG ***
figure(99)

ids1=find(trajs==1);
ids2=find(trajs==2);
ids3=find(trajs==3);
ids4=find(trajs==4);

mu_meas = [];
mu_err = [];
mu_evc = [];
mu_pel = [];
mu_pro = [];

for i =1:4
    ids = find(trajs==i);
    mu_meas = [mu_meas,mean(measures(ids))];
    mu_err = [mu_err,mean(errs(ids))];
    mu_evc = [mu_evc,mean(evcs(ids))];
    mu_pel = [mu_pel,mean(pels(ids))];
    mu_pro = [mu_pro,mean(pros(ids))];
end

plot(mu_meas);
hold on;
plot(mu_err);
plot(mu_evc);
plot(mu_pel);
plot(mu_pro);
hold off;





%% ----------------------------------------
%% plot all the datapoints

%plot ACC response to error
figure(1)
set(gcf,'color','w');
scatter(errs,measures,10,'MarkerFaceColor', ...
        proj.param.plot.white,'MarkerEdgeColor', ...
        proj.param.plot.very_light_grey);
hold on;
y_hat_err = FE.Estimate(1) + FE.Estimate(3)*vseq;
plot(vseq,y_hat_err,'r-','LineWidth',3);

xlim([-1,xmax]);
ylim([ymin,ymax]);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

xlabel('Estimated Error');
ylabel('ACC activation');



%plot ACC response to PRO
figure(2)
set(gcf,'color','w');
scatter(pros,measures,10,'MarkerFaceColor', ...
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

xlabel('Estimated PRO');
ylabel('ACC activation');


%plot ACC response to EVC
figure(3)
set(gcf,'color','w');
scatter(evcs,measures,10,'MarkerFaceColor', ...
        proj.param.plot.white,'MarkerEdgeColor', ...
        proj.param.plot.very_light_grey);
hold on;
y_hat_evc = FE.Estimate(1) + FE.Estimate(7)*vseq;
plot(vseq,y_hat_evc,'r-','LineWidth',3);

xlim([xmin,xmax]);
ylim([ymin,ymax]);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

xlabel('Estimated Q-Value');
ylabel('ACC activation');
