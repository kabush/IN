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
logger(['Comparing IN Cog Control Model Predictions       '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.ctrl.in_acc_activ]);
    eval(['! rm -rf ',proj.path.ctrl.in_acc_activ]);
    disp(['Creating ',proj.path.ctrl.in_acc_activ]);
    eval(['! mkdir ',proj.path.ctrl.in_acc_activ]);
end

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
subjects = [];
% trajs = [];

% debug
accs = [];
dbg_errs = [];

% fits
b0 = [];
b1 = [];


%% ----------------------------------------
%% Transform beta-series into affect series {v,a}
subj_cnt = 0;
for i = 1:1 %numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    %% log processing of subject
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
        gm_nii = load_nii([proj.path.mri.gm_mask,'group_gm_mask.nii']);
        mask = double(gm_nii.img);
        brain_size=size(mask);
        mask = reshape(mask,brain_size(1)*brain_size(2)*brain_size(3),1);
        in_brain=find(mask==1);  

        % Load beta-series
        base_nii = load_untouch_nii([proj.path.betas.fmri_in_beta,subj_study,'_',name,'_lss.nii']);
        brain_size = size(base_nii.img);
        
        % Vectorize the base image
        base_img = vec_img_2d_nii(base_nii);
        base_img = reshape(base_img,brain_size(1)*brain_size(2)*brain_size(3),brain_size(4))';

        %% ----------------------------------------
        %% Find ACC activation trajectory
        base_act = base_img(:,in_brain);
        
        %% ----------------------------------------
        %% Load computational models

        %error
        load([proj.path.ctrl.in_err_mdl,subj_study,'_',name,'_mdls.mat']);
        err_mdls = mdls;

        %% %conflict
        %% load([proj.path.ctrl.in_cnf_mdl,subj_study,'_',name,'_mdls.mat']);
        %% cnf_mdls = mdls;
        %% 
        %% %prediction err likelihood
        %% load([proj.path.ctrl.in_pel_mdl,subj_study,'_',name,'_mdls.mat']);
        %% pel_mdls = mdls;
        %% 
        %% %prediction response outcome
        %% load([proj.path.ctrl.in_pro_mdl,subj_study,'_',name,'_mdls.mat']);
        %% pro_mdls = mdls;
        %% 
        %% %expected value of control
        %% load([proj.path.ctrl.in_evc_mdl,subj_study,'_',name,'_mdls.mat']);
        %% evc_mdls = mdls;

        %indices (using err as template)
        indx = reshape(err_mdls.v_indx',1,prod(size(err_mdls.v_indx)));
        
        % % true measure
        % act = base_act(indx,:);
        % acts = [acts,zscore(act)];

        % gm activations
        act = base_act(indx,:);

        % true predictors
        err = reshape(sqrt((err_mdls.v_dcmp').^2),1,prod(size(err_mdls.v_dcmp)))';

        % true trajectory
        traj_box = repmat([1:4],30,1);
        traj = reshape(traj_box',1,prod(size(traj_box)))';

        betas = zeros(1,size(act,2));
        for j=1:size(act,2)
            if(mod(j,1000)==0)
               j
            end
            [b,stat] = robustfit(act(:,j),traj);
            b0 = [b0,b(1)];
            b1 = [b1,b(2)];

        end


    end
  
end




% %% ----------------------------------------
% %% Assemble Measures
% measures = double(measures);
% errs = double(errs-mean(errs));
% cnfs = double(cnfs-mean(cnfs));
% pels = double(pels-mean(pels)); 
% pros = double(pros-mean(pros)); 
% evcs = double(evcs-mean(evcs)); 
% subjects = double(subjects);
% 
% %% Group GLMM fit
% tbl = table(measures,errs,cnfs,pels,pros,evcs,subjects,'VariableNames',{'trg', ...
%                     'err','cnf','pel','pro','evc','subj'});;
% 
% % mdl_fe = fitlme(tbl,['trg ~ 1 + err + cnf + pel + pro + evc']);
% % mdl_re= fitlme(tbl,['trg ~ 1 + err + cnf + pel + pro + evc + ' ...
% %                     '(err|subj) + (cnf|subj) + (pel|subj) + ' ...
% %                     '(pro|subj) + (evc|subj)']);
% 
% 
% mdl_fe = fitlme(tbl,['trg ~ 1 + err']);
% mdl_re= fitlme(tbl,['trg ~ 1 + err + (err|subj)']);
% 
% figure(1)
% mu_meas = [];
% mu_err = [];
% mu_evc = [];
% mu_pel = [];
% mu_pro = [];
% mu_cnf = [];
% for i =1:4
%     ids = find(trajs==i);
%     mu_meas = [mu_meas,mean(measures(ids))];
%     mu_err = [mu_err,mean(errs(ids))];
%     mu_evc = [mu_evc,mean(evcs(ids))];
%     mu_pel = [mu_pel,mean(pels(ids))];
%     mu_pro = [mu_pro,mean(pros(ids))];
%     mu_cnf = [mu_cnf,mean(cnfs(ids))];
% end
% plot(mu_meas,'LineWidth',3);%blue
% hold on;
% plot(mu_err,'LineWidth',2); %orange
% plot(mu_evc,'LineWidth',2); %mustard
% plot(mu_pel,'LineWidth',2); %purple
% plot(mu_pro,'LineWidth',2); %green
% plot(mu_cnf,'LineWidth',2); %cyan
% hold off;
% 
% export_fig 'ACC_vs_evaluations_v.png' -r300
% eval(['! mv ',proj.path.code,'ACC_vs_evaluations_v.png ',proj.path.fig]);
% 
% %%Explore random effects across model types
% disp(' ');
% fe_v_re = compare(mdl_fe,mdl_re);
% 
% mdl = mdl_fe;
% if(fe_v_re.pValue<0.05);
%     disp('Random effects matter');
%     mdl = mdl_re;
% else
%     disp('Random effects DO NOT matter');
% end
% 
% disp(' ');
% 
% %% ----------------------------------------
% %% Examine Main Effect
% [~,~,FE] = fixedEffects(mdl);
% [~,~,RE] = randomEffects(mdl);
% 
% for i =2:numel(FE.pValue)
%     if(FE.pValue(i)<0.05)
%         disp(['Fixed Effect [',num2str(i-1),'] is significant']);
%         disp(['  p=',num2str(FE.pValue(i))]);
%     end
% 
% end
% 
% %% ----------------------------------------
% %% compute effect size
% Rsqr = mdl.Rsquared.Ordinary;
% Fsqr = Rsqr/(1-Rsqr);
% logger(['Overall model fit'],proj.path.logfile);
% logger(['  Rsqr=',num2str(Rsqr)],proj.path.logfile);
% logger(['  Fsqr=',num2str(Fsqr)],proj.path.logfile);
% disp(' ');
% 
% % %%% DEBUGGING HERE
% % tmp_img = base_img(indx,in_brain);
% % tmp_acc = mean(accs,2);
% % 
% % tmp_errs = mean(dbg_errs,2);
% % 
% % Nvox = size(tmp_img,2);
% % 
% % all_b = [];
% % for i = 1:Nvox
% % 
% %     if(mod(i,100)==0)
% %         i
% %     end
% % 
% %     [b,stat] = robustfit(tmp_errs,tmp_img(:,i));
% %     all_b = [all_b,b(2)];
% % 
% % end
