function calc_mvpa_ex_gs_refit_perm(proj,good_sbj_ids,all_raw_acc,all_raw_prm,sbj_prms,affect)

grp_corr_acc = [];
grp_incorr_acc = [];

all_corr_ids = {};
all_incorr_ids = {};

grp_corr_prm = [];
grp_incorr_prm = [];

for i=1:numel(good_sbj_ids)

    id = good_sbj_ids(i);

    %% Test sub-selection based on group performance
    val_ids = setdiff(1:numel(good_sbj_ids),id);
    
    %% Test group performance against prior
    corr_ids = [];
    incorr_ids = [];
    refit_ids = [];

    for j=1:size(all_raw_acc,2)

        %% Compute fraction correct
        acc_set = all_raw_acc(val_ids,j);
        img_frac = mean(acc_set(find(acc_set>=0))); 
        
        %% ----------------------------------------
        %% Binomial test (corr/incorrect)
        [phat ci] = binofit(numel(val_ids)/2,numel(val_ids),0.05);

        warning('off','last');

        %% Store correct/incorrect ids
        if(img_frac>ci(2))
            corr_ids = [corr_ids,j];
        end
        if(img_frac<ci(1))
            incorr_ids = [incorr_ids,j];
        end
        
    end
    
    %% Store all ids
    all_corr_ids{i} = corr_ids;
    all_incorr_ids{i} = incorr_ids;
    
    %% Compute adjusted subject classification performance
    grp_corr_acc = [grp_corr_acc;mean(all_raw_acc(i,corr_ids))];
    grp_incorr_acc = [grp_incorr_acc;mean(all_raw_acc(i,incorr_ids))];

    %% Compute adjusted subject permutation performance
    grp_corr_prm = [grp_corr_prm;mean(all_raw_prm(i,corr_ids))];
    grp_incorr_prm = [grp_incorr_prm;mean(all_raw_prm(i,incorr_ids))];

    
end    

%% GROUP-level Permutation testing here (paired t-tests)
[h p ci stat] = ttest(mean(all_raw_acc,2),mean(all_raw_prm,2));
logger(['Grp sig., pre-refit, p=',num2str(p),...
        ', CI=[',num2str(ci(1)),', ',num2str(ci(2)),']'],...
        proj.path.logfile);

%% REFIT performance (and confidence intervals)
[h p acc_ci stats] = ttest(grp_corr_acc);
[h p prm_ci stats] = ttest(grp_corr_prm);
logger(['Grp accuracy, pre-refit=', ...
        num2str(mean(mean(all_raw_acc,2))), ', post-refit=', ...
        num2str(mean(grp_corr_acc)), ', CI=[',...
        num2str(acc_ci(1)),', ',num2str(acc_ci(2)),']'],...
        proj.path.logfile);

%% SUBJECT-level Permutation testing and Refit performance
Nsig = 0;
Ntot = 0;
for i=1:numel(good_sbj_ids)

    id = good_sbj_ids(i);
    
    corr_ids = all_corr_ids{i};
    acc = all_raw_acc(i,corr_ids);
    
    if(affect=='v')
        prm = sbj_prms{id}.smp_prm_v(:,corr_ids);
    else
        prm = sbj_prms{id}.smp_prm_a(:,corr_ids);
    end
        

    % Find single subj decision that survive permutation testing
    Nsig_stim = 0;
    for j = 1:numel(corr_ids)
        Nnull = ceil(0.025*size(prm,1));
        Nextreme = numel(find(prm(:,j)>acc(j)));
        if(Nextreme<Nnull)
            Nsig_stim = Nsig_stim + 1;
        end
    end
    
    % Null for binary decision
    ncorr = numel(corr_ids);
    [a b] = binofit(ncorr/2,ncorr,0.05);

    warning('off','last');
    
    
    if(Nsig_stim/ncorr>b(2))
        Nsig = Nsig+1;
    end
    
    Ntot = Ntot + 1;
 
end 

logger(['Sing. Subj significant post-refit=',num2str(Nsig),'/',...
        num2str(Ntot)],proj.path.logfile);


