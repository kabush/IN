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
logger(['*************************************************'], ...
       proj.path.logfile);
logger(['Analyzing ER Skill                               '], ...
       proj.path.logfile);
logger(['*************************************************'], ...
       proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.analysis.dynamics]);
    eval(['! rm -rf ',proj.path.analysis.dynamics]);
    disp(['Creating ',proj.path.analysis.dynamics]);
    eval(['! mkdir ',proj.path.analysis.dynamics]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

figure(1)
set(gcf,'color','w');

all_b_t1 = [];
all_b_t2 = [];

for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    % log analysis of subject
    logger([subj_study,'_',name],proj.path.logfile);

    try

        %% Load IN trajectory structures
        load([proj.path.ctrl.in_ctrl,subj_study,'_',name, ...
              '_prds.mat']);

        if(isfield(prds,'v_dcmp'))

            %% ----------------------------------------
            %% Want to make two tests
            %% 1) does e(t-1) predict e(t)
            %% 2) does [e(t-2),e(t-1)] better predict e(t)
            err = prds.v_dcmp.h-prds.v_dcmp.h(:,1);

            t0_vec = reshape(err(:,6:7)',prod(size(err(:,6:7))),1);
            t1_vec = reshape(err(:,5:6)',prod(size(err(:,5:6))),1);
            t2_vec = reshape(err(:,4:5)',prod(size(err(:,4:5))),1);
            t3_vec = reshape(err(:,3:4)',prod(size(err(:,3:4))),1);
            t4_vec = reshape(err(:,2:3)',prod(size(err(:,2:3))),1);

            %% scatter plot specific points        
            scatter(t2_vec,t1_vec,10,'MarkerFaceColor', ...
                    proj.param.plot.white,'MarkerEdgeColor', ...
                    proj.param.plot.light_grey);
            hold on;

            stats1  = regstats(t0_vec,t1_vec,'linear');
            stats2  = regstats(t0_vec,[t1_vec,t2_vec,t3_vec, ...
                                t4_vec],'linear');

            %%             all_b_t1 = [all_b_t1,
            %%             stats1.tstat_beta];
            %%             all_b_t1 = [all_b_t1,
            %%             stats1.tstat_beta];
            %% 
        else
            disp(['  -Could not find v_dcmp for: ',subj_study,'_', ...
                  name],proj.path.logfile);
        end
        
    catch
        % do nothing
        logger(['  -Could not find/load prds for: ',subj_study,'_', ...
                name],proj.path.logfile);
    end

end


%% ----------------------------------------
%% format figure
xmin = -3;
xmax = 3;
ymin = -3;
ymax = 3;

xlim([xmin,xmax]);
ylim([ymin,ymax]);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

xlabel('err(t-1)');
ylabel('err(t)');
