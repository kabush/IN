%%========================================
%%========================================
%%
%% Keith Bush, PhD (2019)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

function [gamma_opt,frac_opt] = calc_q_param_opt(gamma_set,frac_set,act_err)

%% Find the lowest action error obeying the following constraints:
%% 1) favor low action fraction
%% 2) favor maximal discount (to differentiate from simple error)
%% 3) favor minimal error

%% Find globally minimum error
min_err_gbl = min(min(act_err));

%% Find fraction with min global error
min_err_gbl_not_found = 1;
for i=1:numel(frac_set)
    min_err_frac = min(act_err(:,i));
    if((min_err_frac==min_err_gbl) & min_err_gbl_not_found)
        min_frac_id = i;
        min_err_gbl_not_found = 0;
    end
end

%% Find gamma of column with min global error
min_err_gbl_not_found = 1;
for i=1:numel(gamma_set)

    gamma_id = numel(gamma_set)-(i-1);
    min_err = act_err(gamma_id,min_frac_id);

    if((min_err==min_err_gbl) & min_err_gbl_not_found)
        min_gamma_id = numel(gamma_set)-(i-1);
        min_err_gbl_not_found = 0;
    end
end

%% Extract parameters from matrix ids
gamma_opt = gamma_set(min_gamma_id);
frac_opt = frac_set(min_frac_id);
