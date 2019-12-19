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
logger(['Analyzing VALENCE VR Dynamics                    '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

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

h_trg = [];
h_prd = [];
dh_prd = [];
d2h_prd = [];
sbj_re = [];

for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    % log analysis of subject
    logger([subj_study,'_',name],proj.path.logfile);

    try

        %% Load IN trajectory structures
        load([proj.path.ctrl.in_dyn,subj_study,'_',name, ...
              '_prds.mat']);

        if(isfield(prds,'v_dcmp'))

            %2nd deriv exists for times 3,4,5.
            %% Extract usable times
            h = prds.v_dcmp.h(:,4:6);
            h_n1 = prds.v_dcmp.h(:,3:5);
            dh = prds.v_dcmp.dh(:,2:4); %b/c of pruning for drv
            d2h = prds.v_dcmp.d2h(:,1:3); %b/c of pruning for 2nd drv

            %% Vectorize for modeling
            h_trg = [h_trg;reshape(h',prod(size(h)),1)];
            h_prd = [h_prd;reshape(h_n1',prod(size(h_n1)),1)];
            dh_prd = [dh_prd;reshape(dh',prod(size(dh)),1)];
            d2h_prd = [d2h_prd;reshape(d2h',prod(size(d2h)),1)];
            sbj_re = [sbj_re;repmat(i,prod(size(d2h)),1)];

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

disp('-----------------------------------------------');
disp(' Analyzing dynamics via mixed effects          ');
disp('-----------------------------------------------');
disp(' Predicting x(t) given x(t-1),dx(t-1),d2x(t-1) ');
disp(' ');
disp(' mdl_1 = x(t-1) alone');
disp(' mdl_2 = x(t-1) + dx(t-1)');
disp(' mdl_3 = x(t-1) + dx(t-1) + d2x(t-1)');
disp(' ');
disp(' random effects are modeled subject-wise if valid');
disp(' ');

%% ----------------------------------------
%% Group GLMM fit
trg = double(h_trg);
p1 = double(h_prd-mean(h_prd));
p2 = double(dh_prd-mean(dh_prd));
p3 = double(d2h_prd-mean(d2h_prd));
sub = double(sbj_re);

tbl = table(trg,p1,p2,p3,sub,'VariableNames',{'trg','p1','p2','p3','sub'});

mdl_fe_1 = fitlme(tbl,['trg ~ 1 + p1']);
mdl_fe_2 = fitlme(tbl,['trg ~ 1 + p1 + p2']);
mdl_fe_3 = fitlme(tbl,['trg ~ 1 + p2 + p2 + p3']);

mdl_re_1 = fitlme(tbl,['trg ~ 1 + p1 + (1+p1|sub)']);
mdl_re_2 = fitlme(tbl,['trg ~ 1 + p1 + p2 + (1+p1|sub) + (1+p2|sub)']);
mdl_re_3 = fitlme(tbl,['trg ~ 1 + p1 + p2 + p3 + (1+p1|sub) + (1+p2|sub)+ (1+p3|sub)']);

%%Explore random effects across model types
fe_v_re_1 = compare(mdl_fe_1,mdl_re_1);
fe_v_re_2 = compare(mdl_fe_2,mdl_re_2);
fe_v_re_3 = compare(mdl_fe_3,mdl_re_3);

mdl_1 = mdl_fe_1;
if(fe_v_re_1.pValue<0.05);
    disp('   -mdl_1 random effects matter');
    mdl_1 = mdl_re_1;
end

mdl_2 = mdl_fe_2;
if(fe_v_re_2.pValue<0.05);
    disp('   -mdl_2 random effects matter');
    mdl_2 = mdl_re_2;
end

mdl_3 = mdl_fe_3;
if(fe_v_re_3.pValue<0.05);
    disp('   -mdl_3 random effects matter');
    mdl_3 = mdl_re_3;
end

disp(' ');

%% compare position, position+velocity, position+velocity+accleration
h1v2 = compare(mdl_1,mdl_2);
h2v3 = compare(mdl_2,mdl_3);

mdl = mdl_1;
if(h1v2.pValue<0.05)
    disp('   -dx(t-1) matters');
    mdl = mdl_2;
end

if(h2v3.pValue<0.05)
    disp('   -d2x(t-1) matters');
    mdl = mdl_3;
end

disp(' ');

Rres = sum((mdl.residuals).^2);
Rtot = sum((trg-mean(trg)).^2);

Rsqr = 1-(Rres/Rtot);
Fsqr = Rsqr/(1-Rsqr);

disp(['Rsqr=',num2str(Rsqr)]);
disp(['Fsqr=',num2str(Fsqr)]);

disp(' ');

