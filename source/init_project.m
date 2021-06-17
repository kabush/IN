%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

%% ----------------------------------------
%% Seed random number generator
rng(1,'twister');

%% ----------------------------------------
%% Initialize project param structure
proj = struct;

%% ----------------------------------------
%% Link tools
proj.path.tools.kablab = '/home/kabush/lib/kablab/';
addpath(genpath(proj.path.tools.kablab));

proj.path.tools.scralyze = '/home/kabush/lib/scralyze/';
addpath(genpath(proj.path.tools.scralyze));

proj.path.tools.export_fig = '/home/kabush/lib/export_fig/';
addpath(genpath(proj.path.tools.export_fig));

proj.path.tools.nifti = '/home/kabush/lib/nifti/';
addpath(genpath(proj.path.tools.nifti));

proj.path.tools.approxrl = '/home/kabush/lib/approxrl/';
addpath(genpath(proj.path.tools.approxrl));

proj.path.tools.noboxplot = '/home/kabush/lib/noboxplot/';
addpath(genpath(proj.path.tools.noboxplot));

%% ----------------------------------------
%% Project Flag Definitions
proj.flag.clean_build = 1;

%% ----------------------------------------
%% Project Path Definitions

%% Raw data
proj.path.raw_data = '/raw/bush/';
proj.path.raw_cog = [proj.path.raw_data,'cogbehav/'];
proj.path.raw_physio = 'physio';
proj.path.raw_logs = 'logfiles';
proj.path.raw_tabs = 'tabs';
proj.path.atlas = '/home/kabush/atlas/';
proj.path.demo = 'demo';

%% Workspace
proj.path.home = '/home/kabush/workspace/';
proj.path.name = 'IN';
proj.path.code = [proj.path.home,'code/',proj.path.name,'/'];
proj.path.data = [proj.path.home,'data/',proj.path.name,'/'];
proj.path.log =[proj.path.code,'log/']; 
proj.path.fig = [proj.path.code,'fig/'];

%% Subject Lists
proj.path.subj_list = [proj.path.code,'subj_lists/'];

%% Design path (this is a meta source file)
proj.path.design = [proj.path.code,'design/'];

%% Logging (creates a unique time-stampted logfile)
formatOut = 'yyyy_mm_dd_HH:MM:SS';
t = datetime('now');
ds = datestr(t,formatOut);
proj.path.logfile = [proj.path.log,'logfile_',ds,'.txt'];

%% ----------------------------------------
%% Data Output Directory (All top-level names)
proj.path.system.name = 'system/';
proj.path.analysis.name = 'analysis/';
proj.path.betas.name = 'beta_series/';
proj.path.ctrl.name = 'ctrl/';
proj.path.haufe.name = 'haufe/';
proj.path.mri.name = 'mri/';
proj.path.mvpa.name = 'mvpa/';
proj.path.physio.name = 'physio/';
proj.path.trg.name = 'target/';

%% ----------------------------------------
%% Specific Output Paths

%% System paths
proj.path.sys.subjects = [proj.path.data,proj.path.system.name,'subjects/'];

%% MRI paths
proj.path.mri.mri_clean = [proj.path.data,proj.path.mri.name,'mri_clean/'];
proj.path.mri.gm_mask = [proj.path.data,proj.path.mri.name,'gm_mask/'];

%% Beta-Series Paths
proj.path.betas.fmri_ex_beta = [proj.path.data,proj.path.betas.name,'fmri_ex_beta/'];
proj.path.betas.fmri_in_beta = [proj.path.data,proj.path.betas.name,'fmri_in_beta/'];

%% SCR paths
proj.path.physio.scr_clean = [proj.path.data,proj.path.physio.name,'scr_clean/'];
proj.path.betas.scr_ex_beta = [proj.path.data,proj.path.betas.name,'scr_ex_beta/'];
proj.path.betas.scr_in_beta = [proj.path.data,proj.path.betas.name,'scr_in_beta/'];
proj.path.betas.scr_rest_beta = [proj.path.data,proj.path.betas.name,'scr_rest_beta/'];

%% HR paths 
% Under development within project http://github.com/kabush/HR

%% EMG paths
proj.path.physio.emg_clean = [proj.path.data,proj.path.physio.name,'emg_clean/'];
proj.path.betas.emg_ex_beta = [proj.path.data,proj.path.betas.name,'emg_ex_beta/'];
proj.path.betas.emg_in_beta = [proj.path.data,proj.path.betas.name,'emg_in_beta/'];
proj.path.betas.emg_rest_beta = [proj.path.data,proj.path.betas.name,'emg_rest_beta/'];

%% Target paths
proj.path.trg.ex = [proj.path.data,proj.path.trg.name,'target_ex/'];
proj.path.trg.in = [proj.path.data,proj.path.trg.name,'target_in/'];

%% MVPA paths
proj.path.mvpa.fmri_ex_gs_cls = [proj.path.data,proj.path.mvpa.name,'fmri_ex_gs_cls/'];
proj.path.mvpa.fmri_ex_gm_cls = [proj.path.data,proj.path.mvpa.name,'fmri_ex_gm_cls/'];
proj.path.mvpa.fmri_ex_gm_mdl = [proj.path.data,proj.path.mvpa.name,'fmri_ex_gm_mdl/'];
proj.path.mvpa.fmri_ex_gs_perm = [proj.path.data,proj.path.mvpa.name,'fmri_ex_gs_perm/'];

%% Applying EX models to IN and REST components of the experiment
proj.path.mvpa.fmri_ex_via_ex_gm_mdl = [proj.path.data,proj.path.mvpa.name,'fmri_ex_via_ex_gm_mdl/'];
proj.path.mvpa.fmri_in_via_ex_gm_mdl = [proj.path.data,proj.path.mvpa.name,'fmri_in_via_ex_gm_mdl/'];
proj.path.mvpa.fmri_rest_via_ex_gm_mdl = [proj.path.data,proj.path.mvpa.name,'fmri_rest_via_ex_gm_mdl/'];

%% Measure affect entrainment at rest by physio and fMRI
proj.path.analysis.fmri_rest_entrain = [proj.path.data,proj.path.analysis.name,'fmri_rest_entrain/']
proj.path.analysis.emg_rest_entrain = [proj.path.data,proj.path.analysis.name,'emg_rest_entrain/'];
proj.path.analysis.scr_rest_entrain = [proj.path.data,proj.path.analysis.name,'scr_rest_entrain/'];

%% Secondary replication paths
% secondary replication of Front. in Human Neuro. (2018) paper
proj.path.mvpa.fmri_ex_gs_vs_gm = [proj.path.data,proj.path.mvpa.name,'fmri_ex_gs_vs_gm/'];

% secondary replication of SciReports (2018) paper
proj.path.mvpa.fmri_ex_gm_rgr_scr = [proj.path.data,proj.path.mvpa.name,'fmri_ex_gm_rgr_scr/'];
proj.path.mvpa.fmri_ex_gm_rgr_a = [proj.path.data,proj.path.mvpa.name,'fmri_ex_gm_rgr_a/'];

% secondary replication (not previously reported)
proj.path.mvpa.fmri_ex_gm_rgr_zygo = [proj.path.data,proj.path.mvpa.name,'fmri_ex_gm_rgr_zygo/'];
proj.path.mvpa.fmri_ex_gm_rgr_corr = [proj.path.data,proj.path.mvpa.name,'fmri_ex_gm_rgr_corr/'];

% secondary replication of Psychophysiology (~2020; under review) paper
proj.path.mvpa.fmri_ex_gm_rgr_hr = [proj.path.data,proj.path.mvpa.name,'fmri_ex_gm_rgr_hr/'];
proj.path.mvpa.fmri_ex_gm_rgr_v = [proj.path.data,proj.path.mvpa.name,'fmri_ex_gm_rgr_v/'];

%% Haufe path
proj.path.haufe.fmri_ex_gm_mdl = [proj.path.data,proj.path.haufe.name,'haufe_ex_gm_mdl/'];

%% Intrinsic (IN) control path
proj.path.ctrl.in_dyn = [proj.path.data,proj.path.ctrl.name,'in_dyn/']; ...
proj.path.ctrl.in_ica = [proj.path.data,proj.path.ctrl.name,'in_ica/'];
proj.path.ctrl.in_err_mdl = [proj.path.data,proj.path.ctrl.name,'in_err_mdl/'];
proj.path.ctrl.in_evc_opt_mdl = [proj.path.data,proj.path.ctrl.name,'in_evc_opt_mdl/'];
proj.path.ctrl.in_evc_icv_mdl = [proj.path.data,proj.path.ctrl.name,'in_evc_icv_mdl/'];
proj.path.ctrl.in_pro_opt_mdl = [proj.path.data,proj.path.ctrl.name,'in_pro_opt_mdl/'];

%% Extrinsic (EX) analysis path
proj.path.analysis.gs_cls_refit = [proj.path.data,proj.path.analysis.name,'gs_cls_refit/'];
proj.path.analysis.gs_cls_refit_perm = [proj.path.data,proj.path.analysis.name,'gs_cls_refit_perm/'];
proj.path.analysis.ex_gs_vs_gm = [proj.path.data,proj.path.analysis.name,'ex_gs_vs_gm/'];
proj.path.analysis.ex_scr_a = [proj.path.data,proj.path.analysis.name,'ex_scr_a/'];
proj.path.analysis.ex_emg_v = [proj.path.data,proj.path.analysis.name,'ex_emg_v/'];

%% Intrinsic (IN) analysis path
proj.path.analysis.vr_skill = [proj.path.data,proj.path.analysis.name,'vr_skill/'];
proj.path.analysis.in_cv_cmb_3dlme = [proj.path.data,proj.path.analysis.name,'in_cv_cmb_3dlme/'];
proj.path.analysis.in_cv_cmb_clust_thresh = [proj.path.data,proj.path.analysis.name,'in_cv_cmb_clust_thresh/'];
proj.path.analysis.in_cv_cmb_ccm_effect = [proj.path.data,proj.path.analysis.name,'in_cv_cmb_ccm_effect/'];
proj.path.analysis.in_base_3dlme = [proj.path.data,proj.path.analysis.name,'in_base_3dlme/'];
proj.path.analysis.in_base_clust_thresh = [proj.path.data,proj.path.analysis.name,'in_base_clust_thresh/'];
proj.path.analysis.in_scr = [proj.path.data,proj.path.analysis.name,'in_scr/'];
proj.path.analysis.in_emg = [proj.path.data,proj.path.analysis.name,'in_emg/'];

%% ----------------------------------------
%% Task file nomenclature
proj.path.task.name_id1 = 'Identify_run_1';
proj.path.task.name_id2 = 'Identify_run_2';
proj.path.task.name_rest = 'Rest';

%% ----------------------------------------
%% Project Parameter Definitions

%% Data source
proj.param.studies = {'CTM','INCA'};

%% fMRI Processing param
proj.param.mri.TR = 2.0;
proj.param.mri.slices = 37;
proj.param.mri.slice_pattern = 'seq+z';
proj.param.mri.do_anat = 'yes';
proj.param.mri.do_epi = 'yes';
proj.param.mri.FD_thresh = 0.5;
proj.param.mri.FD_bad_frac = 0.5;

proj.param.mri.tasks = 'identify rest'; % modulate 
proj.param.mri.scans = 'run1 run2';
proj.param.mri.rest_scans = 'run1';

%% *** Annoying extra parameter (silently swear at Philips software
%% engineers) ***  This shift is due to manner in which the design
%% was orginally constructed to accomodate the real-time
%% processing pipeline.  Prior to the Philips R5 upgrade
%% we were dropping 4 inital TRs, so the design built this in.
%% After the R5 upgrade we were dropping zero TRs but the
%% first TR is processed strangely and so is skipped. To
%% adjust for this we shift the design earlier in time by 3*TRs
%% (TR=2s). Basic problem is that the design assumed an 18 s transient period
%% at the start of the identification runs which changed to 12 s
%% following R5 upgrades (shift was introduced to keep original
%% design files intact (possibly bad decision in the long run)
proj.param.trg.r5_shift = -6;

%% Supervised learning labels of stimuli
proj.param.trg.ex_id = 1;
proj.param.trg.in_id = 2;
proj.param.trg.feel_id = 3;

%% values representing binarized valence/arousal classes
proj.param.trg.pos_class = 1;
proj.param.trg.neg_class = -1;

%% Likert scores adjustment parameters
proj.param.trg.dummy_score = -1;
proj.param.trg.mid_score = 5.0; % used to binarize classes

%% Cognitive dynamics labels of stimuli
proj.param.trg.cogdyn.in_id = 1;
proj.param.trg.cogdyn.cue_id = 2;
proj.param.trg.cogdyn.feel_id = 3;
proj.param.trg.cogdyn.rest_id = 4;

%% Cognitive dynamics numbers of stimuli (for analysis)
proj.param.trg.cogdyn.n_stim = 1;
proj.param.trg.cogdyn.n_cue = 1;
proj.param.trg.cogdyn.n_feel = 4;
proj.param.trg.cogdyn.n_rest = 1;

%% Start times of feel TRs relative to IN stimulus times
proj.param.trg.feel_times = 4:proj.param.mri.TR:10;
proj.param.trg.cue_times = 2; 
proj.param.trg.post_in_rest_times = 12;

%% Length of stimulus (in seconds)
proj.param.trg.stim_t = 2;

%% Length of the tasks (in units of TR)
proj.param.mri.n_trs_id1 = 282;
proj.param.mri.n_trs_id2 = 282;
proj.param.mri.n_trs_rest = 225;
proj.param.mri.n_trs_mod1 = 310;
proj.param.mri.n_trs_mod2 = 310;

%% Design construction fidelity (20 hz) 
%% all designs are manufactured at hi-fidelity
%% before downsampling to match fMRI acquisition
%% rate to minimize noise caused by slow TR
proj.param.betas.hirez = 20;

%% Biopac channels
proj.param.physio.chan_hr = 1;
proj.param.physio.chan_rsp = 2;
proj.param.physio.chan_scr = 3;
proj.param.physio.chan_emg_zygo = 4;
proj.param.physio.chan_emg_corr = 5;

%% Biopac recording freqs
proj.param.physio.hz_scr = 2000;
proj.param.physio.hz_emg = 2000;
proj.param.physio.hz_hr = 2000;

%% SCR analysis parameters
proj.param.physio.scr.filt_med_samp = 0.01; %(Bach 2015)
proj.param.physio.scr.filt_high = 0.0159;   %halfway between .05 and
                                            %.0159
                                            %(Staib 2015)
proj.param.physio.scr.filt_low = 5;
proj.param.physio.scr.filt_type = 2; 

%% EMG analysis parameters
proj.param.physio.emg.filt_low = 500.0; %% reference???
proj.param.physio.emg.filt_high = 10.0; %% reference???
proj.param.physio.emg.filt_type = 2;

%% MVPA parameters
proj.param.mvpa.kernel = 'linear';
proj.param.mvpa.n_resamp = 30; 
proj.param.mvpa.n_perm = 1000;

%% REST parameters (for entrainment calcs)
proj.param.rest.n_pseudo = 100;
proj.param.rest.n_resample = 30;
proj.param.rest.n_trs_trans = 5;
proj.param.rest.n_trs_tail = 10;

%% EVC gridsearch parameters
proj.param.ctrl.discount_set = [0:.1:1];
proj.param.ctrl.reward_frac_set = [0:.2:0.4];
% the full frac_set [0:0.2:1.0] was previously run (Jan 27, 2020)
% and no change in parames exists above 0.2.  Running 0-0.4 for
% completeness (both Valence and Arousal)

%% Control analysis variable names
proj.param.ctrl.ccm_z_names = {'aff','traj','err','pro','evc','age','sex','yint',...
                    'aff_sex','traj_sex','err_sex','pro_sex','evc_sex'};
proj.param.ctrl.ccm_z_ids = {25,27,29,31,33,35,37,39,41,43,45,47,49};

proj.param.ctrl.ccm_f_names = {'yint','aff','sex','age','traj','err','pro','evc',...
                    'aff_sex','aff_age','sex_age','sex_traj','age_traj','sex_err',...
                    'age_err','sex_pro','age_pro','sex_evc','age_evc','aff_sex_age',...
                    'traj_sex_age','err_sex_age','pro_sex_age','evc_sex_age'};
proj.param.ctrl.ccm_f_ids = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};

%% Base analysis variable names
proj.param.ctrl.base_z_names = {'aff','traj','yint','age','sex','aff_sex','traj_sex'};
proj.param.ctrl.base_z_ids = {13,15,17,19,21,23,25};

proj.param.ctrl.base_f_names = {'aff','traj','yint','age','sex','aff_sex','traj_sex',...
                   'aff_age','traj_age','sex_age','aff_sex_age','traj_sex_age'};
proj.param.ctrl.base_f_ids = {1,4,0,3,2,5,8,6,9,7,10,11};
proj.param.ctrl.ica_ids = 1:18; %%Ray (2013) ICAs to be used

%% Haufe parameters
proj.param.haufe.npermute = 1000;
proj.param.haufe.chunk = 10;

%% Plotting parameters
proj.param.plot.axisLabelFontSize = 18;
proj.param.plot.circleSize = 10;
proj.param.plot.white = [1,1,1];
proj.param.plot.very_light_grey = [.9,.9,.9];
proj.param.plot.light_grey = [.8,.8,.8];
proj.param.plot.dark_grey = [.6,.6,.6];
proj.param.plot.axis_nudge = 0.1;
proj.param.plot.blue = [0,0,1];
proj.param.plot.red = [1,0,0];

%% ----------------------------------------
%% Processing progress and quality control

% Processing flags
proj.process = struct();

proj.process.mri = 0;
proj.process.mask = 0;
proj.process.scr = 0;
proj.process.hr = 0;
proj.process.emg = 0;

proj.process.beta_mri_ex_id = 0;
proj.process.beta_scr_ex_id = 0;
proj.process.beta_hr_ex_id = 0;
proj.process.beta_emg_ex_id = 0;

proj.process.mvpa_ex_gs_cls = 0;
proj.process.mvpa_ex_gm_cls = 0;
proj.process.mvpa_ex_gm_mdl = 0;

% Quality review flags
proj.check = struct();

proj.check.mri = 0;
proj.check.mask = 0;
proj.check.scr = 0;
proj.check.hr = 0;
proj.check.emg = 0;

proj.check.beta_mri_ex_id = 0;
proj.check.beta_scr_ex_id = 0;
proj.check.beta_hr_ex_id = 0;
proj.check.beta_emg_ex_id = 0;

proj.check.mvpa_ex_gs_cls = 0;
proj.check.mvpa_ex_gm_cls = 0;
proj.check.mvpa_ex_gm_mdl = 0;

% Create quality control data structures
subjs = load_subjs(proj);
for i=1:numel(subjs)

    % Create quality control structure
    subj = struct();

    % Assign subject details
    subj.study = subjs{i}.study;
    subj.name = subjs{i}.name;
    subj.id = subjs{i}.id;
    
    % Set master flag
    subj.ok = 1;

    % Assign subject to project process
    proj.process.subjs{i} = subj;

end

%% ----------------------------------------
%% Write out initialized project structure
save('proj.mat','proj');