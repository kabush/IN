%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

function [fstd,fstim,fe,mdl_r2,mdl_p] =  search_ex_emg_v(proj,var_name)
%% ----------------------------------------
%% load subjs
subjs = proj.process.subjs; 

%% ----------------------------------------
%% Load labels;
label_id = load([proj.path.trg.ex,'stim_ids.txt']);
v_score = load([proj.path.trg.ex,'stim_v_scores.txt']);

%% Adjust for extrinsic presentations
v_score = v_score(find(label_id==proj.param.trg.ex_id));

%% Prepare stimulus set search (bipolarity)
std_v = std(v_score);
v_score = v_score - mean(v_score);

fstd_seq = 0:0.1:1.5;

fstd = [];
fstim = [];
fe = [];
mdl_r2 = [];
mdl_p = [];

for j=1:numel(fstd_seq)

    this_fstd = fstd_seq(j);
    keep_ids = find(abs(v_score)>this_fstd*std_v);

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
            load([proj.path.betas.emg_ex_beta,subj_study,'_',name,'_ex_betas.mat']);
        catch
            logger('    Could not find emg beta file for processing.', ...
                   proj.path.logfile);
        end
        
        emg_betas = [eval(['ex_betas.',var_name,'.id1']);...
                     eval(['ex_betas.',var_name,'.id2'])];
        
        
        emg_betas = zscore(emg_betas(keep_ids));
        emg_v_score = zscore(v_score(keep_ids));

        if(~isempty(emg_betas))
            
            measures = [measures;emg_v_score];
            predictors = [predictors;zscore(emg_betas)];
            subjects = [subjects;repmat(i,numel(emg_betas),1)];
                        
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
        mdl = mdl_re;
    else
        % nothing
    end
    
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
    

    fstd = [fstd,this_fstd];
    fstim = [fstim,numel(keep_ids)/numel(v_score)];
    fe = [fe,FE.Estimate(2)];
    mdl_r2 = [mdl_r2,Rsqr];
    mdl_p = [mdl_p,FE.pValue(2)];
    
end
