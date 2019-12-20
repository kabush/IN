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

% %% ============================================================
% %% PHASE 0: Project Initialization and Preprocessing 
% %% ============================================================
% 
%% ----------------------------------------
%% STEP 1: Initialize the projects directories and parameters.
init_project;

% %%  ----------------------------------------
% %% STEP 2: Clear and reconstruct the project data folder
% clean_project;

% %%  ----------------------------------------
% %% STEP 3: Check demographics balance
% check_demo;

% %%  ----------------------------------------
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
% %%  ----------------------------------------
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
% % - project_summary;   %% master summary 
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
% % % %% ********* VERY SLOW BELOW **********
% % %% ----------------------------------------
% % %% STEP 9: Hyperplane encoding analysis (see SciReports 2018 paper)
% % 
% % % Encodings of Haufe-transformed hyperplanes
% % haufe_fmri_ex_gm; (uses global permuation test ... slow)
% 
% % %% ----------------------------------------
% % %% STEP 10: V vs A hyperplane cosine sim (see SciReports 2018 paper)
% % % (TBD)
% % %% ********* VERY SLOW ABOVE  **********
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
% calc_scr_in_beta; % (pilot)
% calc_emg_in_beta; % (pilot)
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

% %% ============================================================
% %% PHASE 3: Characterizing mFC function
% %% ============================================================

% %% Construct Basic cognitive control models (CCMs)  
% ccm_err_fmri_in_gm;     % error model
% ccm_cnf_fmri_in_gm;   % conflict model (neu vs extrem)
% % ccm_cnf_alt_fmri_in_gm; % conflict model (val vs aro; not using)
% ccm_pel_fmri_in_gm;     % prediction error likelihood
% ccm_pro_fmri_in_gm;     % predicted response outcome

% %% Construct Reinforcement Learning (i.e. EVC) CCM.  State-space
% % is constructed from Ray, 2013 (emotion ICs, 5 of 20) in which
% % mFC has been excluded (mFC will be the CC space).
% reshape_ica;  % ***TICKET*** finalize dACC and IC interaction and ICs

% % Conduct grid search of EVC parm space (mix of err/action-cost)
% ccm_evc_fmri_in_gm_gridsearch; % Q-func. param gridsearch (CNS
%                                % 2020 ***VERY SLOW***)
%                                % Need to select best for next
%                                % step *** TICKET ***
%
% analyze_evc_fmri_in_gm;  % Find best meta-parameters of RL
%
% ccm_evc_fmri_in_gm;  % compute EVC cog mdls w/ max
%                      % params (fit to valence)
% 
% %% Compute CCMs Activations
% analyze_in_fmri_3dlme;
%
% %% Estimate & Apply Cluster Thresholds
% calc_in_clust_thresh_3dlme;
% 
% %% Compute Prediction Effects for CCMs
analyze_ccm_effect;

% %% Compare Predictions Effects
% % (TBD) 
% 
% 
% %% ----------------------------------------
% %% OLD CODE BELOW
% % Predict CCMs from icaACC masked beta-series
% analyze_ccm_ALL_v;
% analyze_ccm_icaACC_traj;
% analyze_ccm_icaACC_v;
% analyze_ccm_icaACC_a;
% 
% % Predict CCMs from icalPFC masked beta-series
% analyze_ccm_icalPFC_v;
% 
% %% ============================================================
% %% PHASE 4: Secondary Validation of IN
% %% ============================================================
% 
% %% ***TICKET*** path to dynamics have changes from data/in_ctrl to
% %% data/in_dyn ... will likely break analysis code below.
% 
% %% ------------------------------------------------------------ 
% %% STEP ???: Analyze IN SCR Response (pilot)
% analyze_in_scr;
% analyze_in_emg;
%
% summarize_vr; % -> this information combined with redcap cogbehav data

toc
