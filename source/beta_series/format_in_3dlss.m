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
    disp(['Removing ',proj.path.trg.in]);
    eval(['! rm -rf ',proj.path.trg.in]);
    disp(['Creating ',proj.path.trg.in]);
    eval(['! mkdir ',proj.path.trg.in]);
end

%% Load designs
load([proj.path.design,'run1_design.mat']);
load([proj.path.design,'run2_design.mat']);

%% Extract Intrinsic Stimulation Times (shifted for R5 upgrade)
run1_in_stim_times = run1_design.in_time_seq'+proj.param.trg.r5_shift;
run2_in_stim_times = run2_design.in_time_seq'+proj.param.trg.r5_shift;

%% Extract Cue Stimulation Times
run1_cue_stim_times = [];
for i=1:numel(run1_in_stim_times)
    run1_cue_stim_times = [run1_cue_stim_times; proj.param.trg.cue_times+run1_in_stim_times(i)];
end

run2_cue_stim_times = [];
for i=1:numel(run1_in_stim_times)
    run2_cue_stim_times = [run2_cue_stim_times; proj.param.trg.cue_times+run2_in_stim_times(i)];
end

%% Extract Feel Stimuluation Times
run1_feel_stim_times = [];
for i=1:numel(run1_in_stim_times)
    run1_feel_stim_times = [run1_feel_stim_times, proj.param.trg.feel_times+run1_in_stim_times(i)];
end

run2_feel_stim_times = [];
for i=1:numel(run2_in_stim_times)
    run2_feel_stim_times = [run2_feel_stim_times, proj.param.trg.feel_times+run2_in_stim_times(i)];
end

%% Extract Rest Times (1 rest immediately following Feel)
run1_rest_stim_times = [];
for i=1:numel(run1_in_stim_times)
    run1_rest_stim_times = [run1_rest_stim_times; proj.param.trg.post_in_rest_times+run1_in_stim_times(i)];
end

run2_rest_stim_times = [];
for i=1:numel(run2_in_stim_times)
    run2_rest_stim_times = [run2_rest_stim_times; proj.param.trg.post_in_rest_times+run2_in_stim_times(i)];
end

%% Build Run 1 TARGETS
run1_stim_times = [run1_in_stim_times;
                   run1_cue_stim_times;
                   run1_feel_stim_times';
                   run1_rest_stim_times];

run1_stim_ids = [repmat(1,numel(run1_in_stim_times),1);
                 repmat(2,numel(run1_cue_stim_times),1);
                 repmat(3,numel(run1_feel_stim_times),1);
                 repmat(4,numel(run1_rest_stim_times),1)];

[y,indices] = sort(run1_stim_times);

s_run1_stim_times = run1_stim_times(indices);
s_run1_stim_ids = run1_stim_ids(indices);

%% Build Run 2 TARGETS
run2_stim_times = [run2_in_stim_times;
                   run2_cue_stim_times;
                   run2_feel_stim_times';
                   run2_rest_stim_times];

run2_stim_ids = [repmat(1,numel(run2_in_stim_times),1);
                 repmat(2,numel(run2_cue_stim_times),1);
                 repmat(3,numel(run2_feel_stim_times),1);
                 repmat(4,numel(run2_rest_stim_times),1)];

[y,indices] = sort(run2_stim_times);

s_run2_stim_times = run2_stim_times(indices);
s_run2_stim_ids = run2_stim_ids(indices);

%% Combine Run 1 & Run2
offset = proj.param.mri.n_trs_id1*proj.param.mri.TR
stim_times = [s_run1_stim_times; s_run2_stim_times+offset];
stim_ids = [s_run1_stim_ids;s_run2_stim_ids];

%% Write out targets and labels
filename = [proj.path.trg.in,'stim_times.1D'];
fid = fopen(filename,'w');
fprintf(fid,'%6.2f\n',stim_times);
fclose(fid);

%% Save out
save([proj.path.trg.in,'stim_ids.txt'],'stim_ids','-ascii')
