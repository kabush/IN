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

% %%  ----------------------------------------
% %% STEP 2: Clear and reconstruct the project data folder
% clean_project;
% 
% %%  ----------------------------------------
% %% STEP 3: Check demographics balance
% check_demo;
% 
% %% ----------------------------------------
% %% STEP 3: Preprocess raw data (wrangling, filtering, formatting)
% 
% %% fMRI data
% preprocess_fmri;
% preprocess_mask;
% 
% %% Physio data
% preprocess_scr; 
% preprocess_emg; % (pilot)
% % - preprocess_hr (TBD); 
% 
% %% Cognitive data
% % - preprocess_cog (TBD);
% 
% %% ----------------------------------------
% %% STEP 4: Run quality check (system) on preprocessing outcomes
% check_mri;
% check_scr;
% % - check_hr; (TBD)
% check_emg;
% % - check_cog (TBD)
% 
% %% ============================================================
% %% PHASE 1: Modeling Affective Brain States
% %% ============================================================
% 
% %% ----------------------------------------
% %% STEP 1: Format Extrinsic Stimuli Design
% format_ex_3dlss; 
% 
% %% ----------------------------------------
% %% STEP 2: Calculate Extrinsic (EX) Stimuli Beta-Series
% 
% %% fMRI data
% calc_fmri_ex_beta;
% 
% %% Physio data
% calc_scr_ex_beta;
% % - calc_hr_ex_beta (TBD);
% % - calc_emg_ex_beta (TBD);
% 
% %% ----------------------------------------
% %% STEP 3: Run quality check (system) on ex_beta series
% check_mri_ex_beta;
% check_scr_ex_beta;
% % - check_hr_ex_beta (TBD);
% % - check_emg_ex_beta (TBD);
% % - project_summary; %% master summary 
% 
% %% ----------------------------------------
% %% STEP 4: Conduct MVPA for Extrinsic Stimuli of Sys. I.D.
% 
% mvpa_fmri_ex_gs_cls; % intra-subj Gram-Schmidt MVPA classification
%                      % performs performance estimation using LOOCV
%                      % basis for stimulus refitting (see below)
% 
% mvpa_fmri_ex_gm_cls; % intra-subj whole-brain GM MVPA classifications
%                      % performs performance estimation using LOOCV
%                      % also constructs and saves single model for
%                      % application to IN formats% 
% 
% mvpa_fmri_ex_gm_mdl; % intra-subj whole-brain GM MVPA models (all data)
% 
% %% ----------------------------------------
% %% STEP 5: Data-driven Analysis of Classification (see Frontiers 2018 paper)
% analyze_mvpa_fmri_ex_gs_cls_refit;
% 
% %% ----------------------------------------
% %% STEP 6: Run quality check (system) on mpva
% check_mvpa_ex_gs_cls;
% check_mvpa_ex_gm_cls;
% check_mvpa_ex_gm_mdl;
% 
% %% ----------------------------------------
% %% STEP 7: compare GS vs GM features (see Frontiers 2018 supplemental)
% analyze_mvpa_fmri_ex_gs_vs_gm;
% 
% %% ----------------------------------------
% %% STEP 8: Analyze EX Physiology Response (compare to brain state)
% 
% % Analyze SCR (see SciReports 2018 paper)
% analyze_ex_gm_scr_a;     % predicting arousal from scr 
% mvpa_fmri_ex_gm_rgr_scr; % predicting scr from brain state
% mvpa_fmri_ex_gm_rgr_a;   % predicting aro from brain state (regr)
% 
% %% % Analyze HR deceleration (see Kayla et al. [2019] paper, in review)
% %% analyze_ex_gm_hr_v; (TBD)
% %% mvpa_fmri_ex_gm_rgr_hr; (TBD)
% %% mvpa_fmri_ex_gm_rgr_v; (TBD)
% 
% % %% ************************************
% % %% ************************************
% % %% ********* VERY SLOW BELOW **********
% % 
% % %% ----------------------------------------
% % %% STEP 9: Hyperplane encoding analysis (see SciReports 2018 paper)
% % 
% % % Encodings of Haufe-transformed hyperplanes
% % haufe_fmri_ex_gm; (uses global permuation test ... slow)
% % 
% % %% ----------------------------------------
% % %% STEP 10: V vs A hyperplane cosine sim (see SciReports 2018 paper)
% % % (TBD)
% % 
% % %% ********* VERY SLOW ABOVE  **********
% % %% *************************************
% % %% *************************************
% 
% %% ============================================================
% %% PHASE 2: Modeling Intrinsic Neuromodulation (IN) of Affect
% %% ============================================================
% 
% %% ----------------------------------------
% %% STEP 1: Format project design for IN afni-based beta-series
% format_in_3dlss;
% 
% %% ----------------------------------------
% %% STEP 2: Calcuate Intrinsic (IN) Stimuli Beta-Series
% 
% %% fMRI data
% calc_fmri_in_beta;
% 
% %% Physio data
% calc_scr_in_beta;
% calc_emg_in_beta;
% % - calc_hr_in_beta (TBD);
% 
% %% ----------------------------------------
% %% STEP 3: Run quality check (system) on beta-series
% % (TBD)
% 
% %% ----------------------------------------
% %% STEP 4: Compute Intrinsic Neuromodulation Dynamics
% dynamics_fmri_in_gm;
% 
% %% ----------------------------------------
% %% STEP 5: Run quality check (system) on dynamics
% % (TBD)
% 
% %% ----------------------------------------
% %% STEP 6: Test significance of VR effect (output figs)
% analyze_in_skill;
% 
% %% ============================================================
% %% PHASE 3: Characterizing mFC function
% %% ============================================================
% 
% %% ------------------------------------------------------------
% %% STEP 1: Construct cognitive control models (CCMs)  
% ccm_err_fmri_in_gm;     % error model
% 
% %% ************************************
% %% ********* VERY SLOW BELOW **********
% 
% % Construct Reinforcement Learning (i.e. EVC) CCM.  State-space
% % is constructed from Ray, 2013 (emotion ICs, 5 of 20) in which
% % mFC has been excluded (mFC will be the CC space).
% reshape_ica;
% 
% % Conduct grid search of EVC parm space (mix of err/action-cost)
% ccm_evc_fmri_in_gm_gridsearch; % Q-func. param gridsearch (CNS
%                                % 2020 ***VERY SLOW***)
%                                % Need to select best for next
%                                % step *** TICKET ***
%                                
% %% ********* VERY SLOW ABOVE **********
% %% ************************************
% 
% %% Select optimal parameters (separately for V & A)w
% select_opt_evc_fmri_in_params;  % Find best meta-parameters of RL
% 
% %% Compute Q-values (separately for V & A)
% ccm_evc_cv_fmri_in_gm; % compute EVC (cross-validated inter-subj)
% 
% %% Compute PRO-values
% ccm_pro_opt_cv_fmri_in_gm;
% 
% %% ------------------------------------------------------------
% %% STEP 2: Identify (and compare) CCM neural correlates
% 
% %% Compute 3dLME Effects (full model with ccms)
% analyze_in_cv_cmb_fmri_3dlme;
% 
% %% Estimate & apply cluster thresholds
% analyze_in_cv_cmb_clust_thresh_3dlme;
% 
% %% Estimate CCM Effects
% analyze_in_cv_cmb_ccm_effect;
% 
% %% ------------------------------------------------------------
% %% STEP 3: Identify BASE TASK neural correlates
% 
% %%Compute 3dLME Effects (base with only trajectory)
% analyze_in_base_fmri_3dlme;

%% Estimate & apply cluster thresholds
analyze_in_base_clust_thresh_3dlme;

% %% ============================================================
% %% PHASE 4: Secondary Validation of IN
% %% ============================================================
% 
% %% ------------------------------------------------------------ 
% %% STEP 1: Validation of volitional affect induction
% analyze_in_scr;
% analyze_in_emg;

toc
