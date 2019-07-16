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
%% ----------------------------------------
%% STEP 3: Run quality check (system) on ex_beta series
check_mri_ex_beta;
check_scr_ex_beta;
% % - check_hr_ex_beta (TBD);
% % - check_emg_ex_beta (TBD);
% % - project_summary;   %% master summary of 
% 
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
%% STEP 6: Run quality check (system) on mpva
check_mvpa_ex_gs_cls;
check_mvpa_ex_gm_cls;
check_mvpa_ex_gm_mdl;

% %% ----------------------------------------
% %% STEP 7: compare GS vs GM features (see Frontiers 2018 supplemental)
% % (TBD)
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
% %% ----------------------------------------
% %% STEP 9: Hyperplane encoding analysis (see SciReports 2018 paper)
% 
% % Encodings of Haufe-transformed hyperplanes
% haufe_fmri_ex_gm; (uses global permuation test ... slow)
% 
% %% ----------------------------------------
% %% STEP 10: V vs A hyperplane cosine sim (see SciReports 2018 paper)
% % (TBD)
% 
% %% ============================================================
% %% PHASE 2: Modeling Intrinsic Neuromodulation of Affect
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
% %% STEP 6: Critically test ACC function (AIM 1)
% 
% % format Ray, 2013 70 ICA component (icaACC) to align with beta-series
% % format Ray, 2013 20 ICA component (RL state) to align with beta-series
% reshape_ica;
% 
% % Compute cognitive control models (CCM) of ACC function
% ccm_err_fmri_in_gm;  %% error model
% ccm_cnf_fmri_in_gm;  %% conflict model
% ccm_pel_fmri_in_gm;  %% prediction error likelihood
% ccm_pro_fmri_in_gm;  %% predicted response outcome
% ccm_evc_fmri_in_gm;    %% expected value of control (Q-value) valence
% analyze_Q_in_vr;       %% ***TICKET*** temporarily here...move down
% 
% % Predict CCMs from icaACC masked beta-series
% analyze_ccm_icaACC_v;
% analyze_ccm_icaACC_a;

% Compare prediction performance
% TBD ((QUESTION: Do we first want to exclude non-performers (Using VR
% Skill via single subj significance?)








% %% ------------------------------------------------------------ 
% %% ------------------------------------------------------------ 
% %%  Steps below are in-progress analysis (pilot data/posters)
% %% ------------------------------------------------------------ 
% %% ------------------------------------------------------------ 

%% ***TICKET*** path to dynamics have changes from data/in_ctrl to
%% data/in_dyn ... will likely break analysis code below.

% %% ------------------------------------------------------------ 
% %% STEP ???: Analyze IN SCR Response (pilot)
% analyze_in_scr;
% analyze_in_emg;

% %% ------------------------------------------------------------ 
% %% Analyze IN VR Cognitive Dynamics
% %% These scripts used GLMM to determine the significance affect
% %% dynamics significantly contribute to predicting future affect
% %% 
% %% TICKET: Current code is clunky, repeating scripts separately for
% %% V and A.
% analyze_v_in_vr_dynamics;
% analyze_a_in_vr_dynamics;
% 
% %% ------------------------------------------------------------ 
% %% Analyze IN VR Control Performance
% %% These scripts determine whether VR succeeded groupwise and which
% %% specific subjects significantly conducted VR.  Then, for the 
% %% significant subjects only, mean control trajectores were
% %% analyzed to understand potential control biases.  Subjects are
% %% labeled if they (on average) lose control in either a positive
% %% or negative affect, where the goal is zero deviation from the
% %% affect induced by the cue image.  Analysis is performed
% %% separately for valence and arousal, respectively. 
% 
% %% TICKET: Current code is clunky, repeating scripts separately for V and A.
% analyze_v_in_vr_skill;
% analyze_v_in_vr_skill_sex_diffs; % (pilot)
% analyze_v_in_vr_labels;
% analyze_v_in_vr_extr_subjs;
% 
% analyze_a_in_vr_skill;
% analyze_a_in_vr_labels;
% analyze_a_in_vr_extr_subjs;
%
% summarize_vr; % -> this information combined with redcap cogbehav data
% 
% %% ------------------------------------------------------------ 
% %% Conduct MVPA for Cognitive Dynamics
% % TBD
% 

% %% IN Dynamics (global permutation test)
% TBD

toc
