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

%% Load the designs
load([proj.path.design,'run1_design.mat']);
load([proj.path.design,'run2_design.mat']);

%% Gather IN normative affect scores
ex_valence = [run1_design.ex_valence_seq,run2_design.ex_valence_seq]
ex_arousal = [run1_design.ex_arousal_seq,run2_design.ex_arousal_seq]

in_valence = [run1_design.in_valence_seq,run2_design.in_valence_seq]
in_arousal = [run1_design.in_arousal_seq,run2_design.in_arousal_seq]


figure(1)
set(gcf,'color','w');

%Underlay extrinsic stimuli (as a reference of the space coverage
scatter(ex_arousal,ex_valence,60,'MarkerFaceColor',...
        proj.param.plot.light_grey,'MarkerEdgeColor',...
        proj.param.plot.light_grey);
hold on;

scatter(in_arousal,in_valence,60,'MarkerFaceColor',...
        proj.param.plot.red,'MarkerEdgeColor',...
        proj.param.plot.red);

xlim([1,8]);
ylim([1,9]);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

export_fig in_stim_affect_summary.png -r300
eval(['! mv ',proj.path.code,'in_stim_affect_summary.png ', ...
      proj.path.fig]);

close all;

%% Gather IN image IDs
iaps_ids = run1_design.img_id_seq(find(run1_design.present_seq==0));
iaps_ids = [iaps_ids,run2_design.img_id_seq(find(run2_design.present_seq==0))];
iaps_ids = sort(iaps_ids);

%% Initialize log section
logger(['*************************************************'],proj.path.logfile);
logger(['IN Stimuli IAPS Image IDs  (Sorted)    '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

for i=1:numel(iaps_ids)
    logger(['  ',num2str(iaps_ids(i))],proj.path.logfile);
end
