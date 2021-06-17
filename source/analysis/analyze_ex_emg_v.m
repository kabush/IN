%%========================================
%%========================================
%%
%% Keith Bush, PhD (2020)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

%% Load in path data
load('proj.mat');

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.analysis.ex_emg_v]);
    eval(['! rm -rf ',proj.path.analysis.ex_emg_v]);
    disp(['Creating ',proj.path.analysis.ex_emg_v]);
    eval(['! mkdir ',proj.path.analysis.ex_emg_v]);
end

%% ----------------------------------------
%% Search through params to test sensitivity of EMG to valence
logger(['*************************************************'],proj.path.logfile);
logger(['Searching EX (Zygomaticus)          '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
[fstd,fstim,z_mdl_fe,z_mdl_r2,z_mdl_p] = search_ex_emg_v(proj,'zygo');

z_emg_search = struct;
z_emg_search.fstd = fstd;
z_emg_search.fstim = fstim;
z_emg_search.z_mdl_fe = z_mdl_fe;
z_emg_search.z_mdl_r2 = z_mdl_r2;
z_emg_search.z_mdl_p = z_mdl_p;
save([proj.path.analysis.ex_emg_v,'z_emg_search.mat'],'z_emg_search');

logger(['*************************************************'],proj.path.logfile);
logger(['Searching EX (Corrugator)          '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
[fstd,fstim,c_mdl_fe,c_mdl_r2,c_mdl_p] = search_ex_emg_v(proj,'corr');

c_emg_search = struct;
c_emg_search.fstd = fstd;
c_emg_search.fstim = fstim;
c_emg_search.c_mdl_fe = c_mdl_fe;
c_emg_search.c_mdl_r2 = c_mdl_r2;
c_emg_search.c_mdl_p = c_mdl_p;
save([proj.path.analysis.ex_emg_v,'c_emg_search.mat'],'c_emg_search');

%% ----------------------------------------
%% Find minimum threshold of  sig. sensitivity
sig_z_ids = find(z_mdl_p<0.05)
sig_c_ids = find(c_mdl_p<0.05)
sig_ids = unique([sig_z_ids,sig_c_ids]);
best_fstd = fstd(sig_ids(1));

%% ----------------------------------------
%% Test EMG measures of valence
logger(['*************************************************'],proj.path.logfile);
logger(['Analyzing EX (Zygomaticus)          '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
calc_ex_emg_v(proj,'zygo',best_fstd);

logger(['*************************************************'],proj.path.logfile);
logger(['Analyzing EX (Corrugator)          '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);
calc_ex_emg_v(proj,'corr',best_fstd);

%% ----------------------------------------
%% Plot summary stats of sensitivity search
figure(1)
set(gcf,'color','w');
plot(fstd,c_mdl_r2,'b-','LineWidth',2)
hold on;
plot(fstd,z_mdl_r2,'r-','LineWidth',2);

csig_id = find(c_mdl_p<0.05);
zsig_id = find(z_mdl_p<0.05);

for i=1:numel(csig_id)
    scatter(fstd(csig_id(i)),c_mdl_r2(csig_id(i)),100,...
            'MarkerFaceColor',proj.param.plot.white,...
            'MarkerEdgeColor',proj.param.plot.blue);
end

for i=1:numel(zsig_id)
    scatter(fstd(zsig_id(i)),z_mdl_r2(zsig_id(i)),100,...
            'MarkerFaceColor',proj.param.plot.white,...
            'MarkerEdgeColor',proj.param.plot.red);
end

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

all_r2 = unique([c_mdl_r2,z_mdl_r2]);
ymax = max(all_r2);
ymin = min(all_r2);
xmax = max(fstd);
xmin = min(fstd);

yrng = ymax-ymin;
xrng = xmax-xmin;
coef = 0.05;

xlim([xmin-coef*xrng,xmax+coef*xrng]);
ylim([ymin-coef*yrng,ymax+coef*yrng]);
box off;

% export hi-def figure
export_fig emg_sensitivity.png -r300
eval(['! mv ',proj.path.code,'emg_sensitivity.png ',proj.path.fig]);

close all;

figure(2)
set(gcf,'color','w');
plot(fstd,fstim,'Color',proj.param.plot.dark_grey,'LineWidth',2)
hold on;

scatter(fstd(sig_ids(1)),fstim(sig_ids(1)),100,...
        'MarkerFaceColor',proj.param.plot.white,...
        'MarkerEdgeColor',proj.param.plot.dark_grey);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

xlim([xmin-coef*xrng,xmax+coef*xrng]);
box off;

% export hi-def figure
export_fig emg_stimset_size.png -r300
eval(['! mv ',proj.path.code,'emg_stimset_size.png ',proj.path.fig]);

close all;
