%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

%% Load in path data
load('proj.mat');

%% Set-up Directory Structure for fMRI Identify run targets
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.trg.ex]);
    eval(['! rm -rf ',proj.path.trg.ex]);
    disp(['Creating ',proj.path.trg.ex]);
    eval(['! mkdir ',proj.path.trg.ex]);
end

%% Load designs
load([proj.path.design,'run1_design.mat']);
load([proj.path.design,'run2_design.mat']);

%% Extract Extrinsic Stimulation Times (shifted for R5 upgrade)
run1_ex_stim_times = run1_design.ex_time_seq'+proj.param.trg.r5_shift; 
run2_ex_stim_times = run2_design.ex_time_seq'+proj.param.trg.r5_shift;

%% Extract Intrinsic Stimulation Times (shifted for R5 upgrade)
run1_in_stim_times = run1_design.in_time_seq'+proj.param.trg.r5_shift;
run2_in_stim_times = run2_design.in_time_seq'+proj.param.trg.r5_shift;

%% Extract Feel Stimuluation Times
run1_feel_stim_times = [];
for i=1:numel(run1_in_stim_times)
    run1_feel_stim_times = [run1_feel_stim_times,proj.param.trg.feel_times+run1_in_stim_times(i)];
end

run2_feel_stim_times = [];
for i=1:numel(run2_in_stim_times)
    run2_feel_stim_times = [run2_feel_stim_times,proj.param.trg.feel_times+run2_in_stim_times(i)];
end

%% ----------------------------------------
%% Build Run 1 TARGETS

% concatenate the times
run1_stim_times = [run1_ex_stim_times;
                   run1_in_stim_times;
                   run1_feel_stim_times'];

% concatenate the ids
run1_stim_ids = [repmat(proj.param.trg.ex_id,numel(run1_ex_stim_times),1);
                 repmat(proj.param.trg.in_id,numel(run1_in_stim_times),1);
                 repmat(proj.param.trg.feel_id,numel(run1_feel_stim_times),1)];

% concatenate the valence scores
run1_stim_v_scores = [run1_design.ex_valence_seq';
                 run1_design.in_valence_seq';
                 repmat(proj.param.trg.dummy_score,numel(run1_feel_stim_times),1)];

% concatenate the arousal scores
run1_stim_a_scores = [run1_design.ex_arousal_seq';
                 run1_design.in_arousal_seq';
                 repmat(proj.param.trg.dummy_score,numel(run1_feel_stim_times),1)];

% assign binary labels to the valence stimuli
run1_stim_v_labs = 0*run1_stim_v_scores; % allocate zero vector
run1_stim_v_labs(find(run1_stim_v_scores>=proj.param.trg.mid_score)) = proj.param.trg.pos_class;
run1_stim_v_labs(find(run1_stim_v_scores<proj.param.trg.mid_score)) = proj.param.trg.neg_class;
run1_stim_v_labs(find(run1_stim_ids==proj.param.trg.feel_id)) = 0; %0 means N/A

% assign binary labels to the arousal stimuli
run1_stim_a_labs = 0*run1_stim_a_scores; % allocate zero vector
run1_stim_a_labs(find(run1_stim_a_scores>=proj.param.trg.mid_score)) = proj.param.trg.pos_class;
run1_stim_a_labs(find(run1_stim_a_scores<proj.param.trg.mid_score)) = proj.param.trg.neg_class;
run1_stim_a_labs(find(run1_stim_ids==proj.param.trg.feel_id)) = 0; %0 means N/A

% order the indices by stimulus time
[y,indices] = sort(run1_stim_times);

% reorder all data by time-sorted indices
s_run1_stim_times = run1_stim_times(indices);
s_run1_stim_ids = run1_stim_ids(indices);
s_run1_stim_v_scores = run1_stim_v_scores(indices);
s_run1_stim_a_scores = run1_stim_a_scores(indices);
s_run1_stim_v_labs = run1_stim_v_labs(indices);
s_run1_stim_a_labs = run1_stim_a_labs(indices);


%% ----------------------------------------
%% Build Run 2 TARGETS (same as Run 1)

% concatenate the times
run2_stim_times = [run2_ex_stim_times;
                   run2_in_stim_times;
                   run2_feel_stim_times'];

% concatenate the ids
run2_stim_ids = [repmat(proj.param.trg.ex_id,numel(run2_ex_stim_times),1);
                 repmat(proj.param.trg.in_id,numel(run2_in_stim_times),1);
                 repmat(proj.param.trg.feel_id,numel(run2_feel_stim_times),1)];

% concatenate the valence scores
run2_stim_v_scores = [run2_design.ex_valence_seq';
                 run2_design.in_valence_seq';
                 repmat(proj.param.trg.dummy_score,numel(run2_feel_stim_times),1)];

% concatenate the arousal scores
run2_stim_a_scores = [run2_design.ex_arousal_seq';
                 run2_design.in_arousal_seq';
                 repmat(proj.param.trg.dummy_score,numel(run2_feel_stim_times),1)];

% assign binary labels to the valence stimuli
run2_stim_v_labs = 0*run2_stim_v_scores; % allocate zero vector
run2_stim_v_labs(find(run2_stim_v_scores>=proj.param.trg.mid_score)) = proj.param.trg.pos_class;
run2_stim_v_labs(find(run2_stim_v_scores<proj.param.trg.mid_score)) = proj.param.trg.neg_class;
run2_stim_v_labs(find(run2_stim_ids==proj.param.trg.feel_id)) = 0; %0 means N/A

% assign binary labels to the arousal stimuli
run2_stim_a_labs = 0*run2_stim_a_scores; % allocate zero vector
run2_stim_a_labs(find(run2_stim_a_scores>=proj.param.trg.mid_score)) = proj.param.trg.pos_class;
run2_stim_a_labs(find(run2_stim_a_scores<proj.param.trg.mid_score)) = proj.param.trg.neg_class;
run2_stim_a_labs(find(run2_stim_ids==proj.param.trg.feel_id)) = 0; %0 means N/A

% order the indices by stimulus time
[y,indices] = sort(run2_stim_times);

% reorder all data by time-sorted indices
s_run2_stim_times = run2_stim_times(indices);
s_run2_stim_ids = run2_stim_ids(indices);
s_run2_stim_v_scores = run2_stim_v_scores(indices);
s_run2_stim_a_scores = run2_stim_a_scores(indices);
s_run2_stim_v_labs = run2_stim_v_labs(indices);
s_run2_stim_a_labs = run2_stim_a_labs(indices);

%% ----------------------------------------
%% Combine Run 1 & Run2
offset = proj.param.mri.n_trs_id1*proj.param.mri.TR; 
stim_times = [s_run1_stim_times; s_run2_stim_times+offset];
stim_ids = [s_run1_stim_ids;s_run2_stim_ids];
stim_v_scores = [s_run1_stim_v_scores;s_run2_stim_v_scores];
stim_a_scores = [s_run1_stim_a_scores;s_run2_stim_a_scores];
stim_v_labs = [s_run1_stim_v_labs;s_run2_stim_v_labs];
stim_a_labs = [s_run1_stim_a_labs;s_run2_stim_a_labs];

%% ----------------------------------------
%% Write out stim times (post-hoc analysis)
filename = [proj.path.trg.ex,'stim_times.1D'];
fid = fopen(filename,'w');
fprintf(fid,'%6.2f\n',stim_times);
fclose(fid);

%% ----------------------------------------
%% Write out the ids, scores and labels
save([proj.path.trg.ex,'stim_ids.txt'],'stim_ids','-ascii');
save([proj.path.trg.ex,'stim_v_scores.txt'],'stim_v_scores','-ascii');
save([proj.path.trg.ex,'stim_a_scores.txt'],'stim_a_scores','-ascii');
save([proj.path.trg.ex,'stim_v_labs.txt'],'stim_v_labs','-ascii');
save([proj.path.trg.ex,'stim_a_labs.txt'],'stim_a_labs','-ascii');
