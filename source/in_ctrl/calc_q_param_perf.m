%%========================================
%%========================================
%%
%% Keith Bush, PhD (2019)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

% Observe Q-function parameters 
function [all_q_perf,all_act_err,all_sig_test] = calc_q_param_perf(proj,...
                                                  discount_set,...
                                                  reward_frac_set,...
                                                  Q_traj_all,...
                                                  Q_rand_all,...
                                                  act_err_all)

%% Calculate set-sizes
[Ndsct,Nfrac,Nsbj,Nvr,Nt] = size(Q_traj_all);

%% META-LOOPS HERE
all_q_traj_dsct = zeros(Ndsct,Nfrac,Nsbj);
all_q_rand_dsct = zeros(Ndsct,Nfrac,Nsbj);
all_q_diff_dsct = zeros(Ndsct,Nfrac,Nsbj);
all_q_perf_dsct = zeros(Ndsct,Nfrac,Nsbj);
act_err_dsct = zeros(Ndsct,Nfrac,Nsbj);

for a=1:Ndsct
    
    for b=1:Nfrac
            
        gamma = discount_set(a);
        rwrd_f = reward_frac_set(b);
        
        all_q_diff = [];
        all_q_rand = [];
        all_q_effx = [];
        
        for c = 1:Nsbj

            %% Q analysis
            q_traj = squeeze(Q_traj_all(a,b,c,:,:));
            q_rand = squeeze(Q_rand_all(a,b,c,:,:));
            
            q_diff = (q_traj-q_rand);
            q_perf = q_diff./abs(q_rand);
            q_diff = reshape(q_diff,prod(size(q_diff)),1);
            q_effx = q_diff/std(vec(q_diff));
            
            all_q_traj_dsct(a,b,c) = median(median(q_traj));
            all_q_rand_dsct(a,b,c) = median(median(q_rand));
            all_q_diff_dsct(a,b,c) = median(median(q_diff));
            all_q_perf_dsct(a,b,c) = median(median(q_perf));
            
            %% Action analysis
            act_err_dsct(a,b,c) = mean(mean(abs(squeeze(act_err_all(a,b,c,:,:)))));
            
        end
        
    end
end

%% observe best parameter combination
all_q_perf = zeros(Ndsct,Nfrac);
all_act_err = zeros(Ndsct,Nfrac);
all_sig_test = zeros(Ndsct,Nfrac);

for a=1:Ndsct
    for b=1:Nfrac
        all_q_perf(a,b) = median(squeeze(all_q_perf_dsct(a,b,:)));
        all_act_err(a,b) = mean(squeeze(act_err_dsct(a,b,:)));
        all_sig_test(a,b) = signrank(squeeze(all_q_perf_dsct(a,b,:)));
    end
end
