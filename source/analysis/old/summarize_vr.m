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

%% Initialize log section
logger(['*************************************************'],proj.path.logfile);
logger([' Organize VR skills for Cog-Behav analysis       '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);


%% ----------------------------------------
%% load analysis structures (VALENCE)
load([proj.path.analysis.vr_skill,'v_sig_subjs.mat']);
load([proj.path.analysis.vr_skill,'v_non_subjs.mat']);
v_sig_subjs = sig_subjs;
v_non_subjs = non_subjs;

%% ----------------------------------------
%% write out analysis summary
fid = fopen([proj.path.analysis.vr_skill,'v_summary.txt'],'w');
fprintf(fid,'%6s %6s %6s %6s\n','Study','Subj','Sig?','Label');
for i=1:numel(v_sig_subjs)
    subj = v_sig_subjs{i};
    fprintf(fid,'%6s %6s %6s %6s\n',subj.study,subj.name,num2str(1),num2str(subj.in_ctrl));
end

for i=1:numel(v_non_subjs)
    subj = v_non_subjs{i};
    fprintf(fid,'%6s %6s %6s %6s\n',subj.study,subj.name,num2str(0),'n/a');
end
fclose(fid);

%% ----------------------------------------
%% load analysis structures (AROUSAL)

load([proj.path.analysis.vr_skill,'a_sig_subjs.mat']);
load([proj.path.analysis.vr_skill,'a_non_subjs.mat']);
a_sig_subjs = sig_subjs;
a_non_subjs = non_subjs;

%% ----------------------------------------
%% write out analysis summary
fid = fopen([proj.path.analysis.vr_skill,'a_summary.txt'],'w');
fprintf(fid,'%6s %6s %6s %6s\n','Study','Subj','Sig?','Label');
for i=1:numel(a_sig_subjs)
    subj = a_sig_subjs{i};
    fprintf(fid,'%6s %6s %6s %6s\n',subj.study,subj.name,num2str(1),num2str(subj.in_ctrl));
end

for i=1:numel(a_non_subjs)
    subj = a_non_subjs{i};
    fprintf(fid,'%6s %6s %6s %6s\n',subj.study,subj.name,num2str(0),'n/a');
end
fclose(fid);


