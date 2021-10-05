%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

tic

%% ----------------------------------------
%% Clean up matlab environment
matlab_reset;

%% ----------------------------------------
%% Link all source code
addpath(genpath('./source/'));

%% ============================================================
%% PHASE 0: Project Initialization and Preprocessing 
%% ============================================================

%% ----------------------------------------
%% STEP 1: Initialize the projects directories and parameters.
init_project;

%%  ----------------------------------------
%% STEP 2: Clear and reconstruct the project data folder
clean_project;

%%  ----------------------------------------
%% STEP 3: Check demographics balance
check_demo;

%% ----------------------------------------
%% STEP 3: Preprocess raw data (wrangling, filtering, formatting)

%% fMRI data
preprocess_fmri;
preprocess_mask;

%% Physio data
preprocess_scr; 
preprocess_emg;

%% ----------------------------------------
%% STEP 4: Run quality check (system) on preprocessing outcomes
check_mri;
check_scr;
check_emg;

%% ============================================================
%% PHASE 1: Modeling Affective Brain States
%% ============================================================

%% ----------------------------------------
%% STEP 1: Format Extrinsic Stimuli Design
format_ex_3dlss; 

%% ----------------------------------------
%% STEP 2: Calculate Extrinsic (EX) Stimuli Beta-Series

%% fMRI data
calc_fmri_ex_beta;

%% Physio data
calc_scr_ex_beta;
calc_emg_ex_beta;

%% ----------------------------------------
%% STEP 3: Run quality check (system) on ex_beta series
check_mri_ex_beta;
check_scr_ex_beta;

%% ----------------------------------------
%% STEP 4: Conduct MVPA for Extrinsic Stimuli of Sys. I.D.

mvpa_fmri_ex_gs_cls; % intra-subj Gram-Schmidt MVPA classification
                     % performs performance estimation using LOOCV
                     % basis for stimulus refitting (see below)

mvpa_fmri_ex_gm_cls; % intra-subj whole-brain GM MVPA classifications
                     % performs performance estimation using LOOCV
                     % also constructs and saves single model for
                     % application to IN formats% 

mvpa_fmri_ex_gm_mdl; % intra-subj whole-brain GM MVPA models (all data)

%% ----------------------------------------
%% STEP 5: Data-driven Analysis of Classification (see Frontiers 2018 paper)
analyze_mvpa_fmri_ex_gs_cls_refit;

%% ----------------------------------------
%% STEP 6: Permutation Analysis of Classification (see biorXiv 2021 revision).
%% Data-driven Analysis of Classification (see Frontiers 2018 paper)
mvpa_fmri_ex_gs_perm; % Permutation testing of GS MVPA classification

%% analysis code here to correct
analyze_mvpa_fmri_ex_gs_cls_refit_perm;

%% ----------------------------------------
%% STEP 7: Run quality check (system) on mpva
check_mvpa_ex_gs_cls;
check_mvpa_ex_gm_cls;
check_mvpa_ex_gm_mdl;

%% ----------------------------------------
%% STEP 8: compare GS vs GM features (see Frontiers 2018 supplemental)
analyze_mvpa_fmri_ex_gs_vs_gm;

%% ----------------------------------------
%% STEP 9: Analyze EX Physiology Response (compare to brain state)

%% Analyze SCR (see SciReports 2018 paper)
analyze_ex_scr_a;        % predicting arousal from scr 
mvpa_fmri_ex_gm_rgr_scr; % predicting scr from brain state
mvpa_fmri_ex_gm_rgr_a;   % predicting aro from brain state (regr)

%% Analyze EMG (see SciReports 2020 paper)
analyze_ex_emg_v;
 
%%  Analyze HR deceleration (see SciReports 2020 paper)
% (see github.com/kabush/HR for code)

%% ************************************
%% ************************************
%% ********* VERY SLOW BELOW **********

%% ----------------------------------------
%% STEP 10: Hyperplane encoding analysis (see SciReports 2018 paper)
haufe_fmri_ex_gm; % (global permuation test null hypothesis)

%% ********* VERY SLOW ABOVE ***********
%% *************************************
%% *************************************

%% ============================================================
%% PHASE 2: Modeling Intrinsic Neuromodulation (IN) of Affect
%% ============================================================

%% ----------------------------------------
%% STEP 1: Format project design for IN afni-based beta-series
format_in_3dlss;

%% ----------------------------------------
%% STEP 2: Calcuate Intrinsic (IN) Stimuli Beta-Series

%% fMRI data
calc_fmri_in_beta;

%% ----------------------------------------
%% STEP 3: Validate EX & IN predictions (Platt scaled)
mvpa_fmri_ex_via_ex_gm_mdl; 
mvpa_fmri_in_via_ex_gm_mdl; 

%% ----------------------------------------
%% STEP 4: Compute Intrinsic Neuromodulation Dynamics
dynamics_fmri_in_gm;

%% ----------------------------------------
%% STEP 5: Test significance of Cued-Recall
analyze_in_skill;

%% ----------------------------------------
%% STEP 6: Predict REST affect entrainment

%% Predict fMRI REST affect
mvpa_fmri_rest_via_ex_gm_mdl;

%% Analyze fMRI REST entrainment (compare to IN)
analyze_fmri_rest_entrain;


%% ============================================================
%% PHASE 3: Characterizing mFC function
%% ============================================================

%% ------------------------------------------------------------
%% STEP 1: Construct cognitive control models (CCMs)  

%% error model (used in evc, cnf & pro)
ccm_err_fmri_in_gm;  

%% Construct Reinforcement Learning (i.e. EVC) CCM.  State-space
%% is constructed from Ray, 2013 (emotion ICs, 18 of 20) in which
%% mFC has been excluded (mFC will be the CC space).
reshape_ica;

%% ************************************
%% ************************************
%% ********* VERY SLOW BELOW **********

%% Conduct grid search of EVC parm space (mix of err/action-cost)
ccm_evc_fmri_in_gm_gridsearch; % Q-func. param gridsearch (CNS
                               % 2020 ***VERY SLOW***)

%% ********* VERY SLOW ABOVE **********
%% ************************************
%% ************************************

%% Select optimal parameters (separately for V & A)
select_opt_evc_fmri_in_params;  % Find best meta-parameters of RL
                                % using CV

%% Compute Q-values (separately for V & A using best params)
ccm_evc_cv_fmri_in_gm; % compute EVC (cross-validated inter-subj)

%% Compute PRO-values
ccm_pro_opt_cv_fmri_in_gm;

%% ------------------------------------------------------------
%% STEP 2: Identify CCM neural correlates

%% Compute 3dLME Effects (full model with ccms)
analyze_in_cv_cmb_fmri_3dlme;

%% Estimate & apply cluster thresholds
analyze_in_cv_cmb_clust_thresh_3dlme;

%% ------------------------------------------------------------
%% STEP 3: Identify BASE TASK neural correlates

%%Compute 3dLME Effects (base with only trajectory)
analyze_in_base_fmri_3dlme;
  
%% Estimate & apply cluster thresholds
analyze_in_base_clust_thresh_3dlme;

%% ============================================================
%% PHASE 4: Secondary Validation of IN
%% ============================================================

%% Calculate IN Physio data
calc_scr_in_beta;
calc_emg_in_beta;

%% Validation of IN induction via physio
analyze_in_scr;
analyze_in_emg;

%% Calculate REST Physio data
calc_scr_rest_beta;
calc_emg_rest_beta;

%% Control for physio entrainment effects (REST vs IN)
analyze_rest_scr_entrain;
analyze_rest_emg_entrain;

%% Display IN norm affect scores, and IAPS IDs
analyze_in_stim;

toc
