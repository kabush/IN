%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

%% ------------------------------------------------------------
%% Clean up matlab environment
matlab_reset;

tic

%% ------------------------------------------------------------
%% Link all source code
addpath(genpath('./source/'));

%% ------------------------------------------------------------
%% STEP 1: Initialize the projects directories and parameters.
init_project;

% %% ------------------------------------------------------------
% %% STEP 2: Clear and reconstruct the project data folder
% clean_project;
% 
% %% ------------------------------------------------------------
% %% STEP 3: Preprocess raw data (wrangling, filtering, formatting)
% 
% %% fMRI data
% preprocess_fmri;
% preprocess_mask;
% % - preprocess_fd (TBD);
% 
% %% Physio data
% preprocess_scr; 
% % - preprocess_hrv (TBD); 
% % - preprocess_emg (TBD);
%  
% %% Cognitive data
% % - preprocess_cog (TBD);
% 
% %% ------------------------------------------------------------
% %% STEP 4: Format Extrinsic Stimuli Design
% format_ex_3dlss; 
% 
% %% ------------------------------------------------------------
% %% STEP 5: Calculate Extrinsic (EX) Stimuli Beta-Series
% 
% %% fMRI data
% calc_fmri_ex_beta;
% 
% %% Physio data
% calc_scr_ex_beta;
% % - calc_hrv_ex_beta (TBD);
% % - calc_emg_ex_beta (TBD);
% 
% %% ------------------------------------------------------------
% %% STEP 6: Conduct MVPA for Extrinsic Stimuli of Sys. I.D.
%  
% %% Classification of Affect Scores
% mvpa_fmri_ex_gs_cls % intra-subj Gram-Schmidt MVPA classification
%                     % performs performance estimation using LOOCV
%                     % basis for stimulus refitting (see below)
% 
% mvpa_fmri_ex_gm_cls % intra-subj whole-brain GM MVPA classifications
%                     % performs performance estimation using LOOCV
%                     % also constructs and saves single model for
%                     % application to IN formats
% 
% mvpa_fmri_ex_gm_mdl % intra-subj whole-brain GM MVPA models
% 
% % ----------------------------------------
% % TICKET  Modify the above mvpa codes to save
% % out the basis function for project to low-dim
% % space.  Will use the inverse to project Haufe-transformed
% % hyperplanes back into GM space for viewing.
% 
% %% ------------------------------------------------------------ 
% %% STEP 7: compare GS vs GM features (Frontiers 2018 paper)
% % TBD
% 
% %% ------------------------------------------------------------ 
% %% STEP 8: Stimulus adjusted performance (see Frontiers 2018 paper)
% % TBD
% 
% %% ------------------------------------------------------------ 
% %% STEP 9: Analyze EX SCR Response (see SciReports 2018 paper)
% analyze_ex_scr; % (draft code, convert to GLMM)
% 
% %% ------------------------------------------------------------ 
% %% STEP 10: V vs A hyperplane cosine sim (see SciReports 2018 paper)
% % TBD
% 
% %% ------------------------------------------------------------
% %% STEP 11: Conduct MVPA for Secondary Measures
% % mvpa_fmri_ex_rgr_scr % (***unworking DRAFT***)
% 
% %% ------------------------------------------------------------ 
% %% STEP 12: Format project design for IN afni-based beta-series
% format_in_3dlss;
% 
% %% ------------------------------------------------------------
% %% STEP 13: Calcuate Intrinsic (IN) Stimuli Beta-Series
% 
% %% fMRI data
% calc_fmri_in_beta;
% 
% %% Physio data
% %- calc_scr_in_beta (TBD);
% 
% %% ------------------------------------------------------------ 
% %% STEP 14: Compute IN VR Cognitive Dynamics
% dynamics_fmri_in_gm;
% 
% %% ------------------------------------------------------------ 
% %% STEP 15: Analyze IN VR Cognitive Dynamics
% %% These scripts used GLMM to determine the significance affect
% %% dynamics significantly contribute to predicting future affect
% %% 
% %% TICKET: Current code is clunky, repeating scripts separately for
% %% V and A.
% analyze_v_in_vr_dynamics;
% analyze_a_in_vr_dynamics;
% 
% %% ------------------------------------------------------------ 
% %% STEP 16: Analyze IN VR Control Performance
% %% These scripts determine whether VR succeeded groupwise and which
% %% specific subjects significantly conducted VR.  Then, for the 
% %% significant subjects only, mean control trajectores were
% %% analyzed to understand potential control biases.  Subjects are
% %% labeled if they (on average) lose control in either a positive
% %% or negative affect, where the goal is zero deviation from the
% %% affect induced by the cue image.  Analysis is performed
% %% separately for valence and arousal, respectively. 
% %%
% %% TICKET: Current code is clunky, repeating scripts separately for V and A.
% analyze_v_in_vr_skill;
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
% %% STEP 17: Conduct MVPA for Cognitive Dynamics
% % TBD
% 
% %% ------------------------------------------------------------ 
% %% STEP 18: Hyperplane analysis (see SciReports 2018 paper)
% %% (Last position b/c it's slow ... lots of resampling)
% 
% %% EX Hyperplanes (global permutation test)
% haufe_fmri_ex_gm;
% 
% %% IN Dynamics (global permutation test)
% TBD

toc
