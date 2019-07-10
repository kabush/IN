%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

%% Load in path data from command line
addpath(home_path);
load('proj.mat');

filein=[proj.path.code,'tmp/var_names'];
fid=fopen(filein,'r');
vars=textscan(fid,'%s');
fclose(fid);

%% Load params of current data to process motion
study=char(vars{1}(1));
subj=char(vars{1}(2));
task=char(vars{1}(3));
scan=char(vars{1}(4));

%% Load motion file
filename = [proj.path.mri.mri_clean,study,'_',subj,'/',task, ...
            '/',scan,'/',study,'.',subj,'.',task,'.',scan,'.motion.1D'];

motion=load(filename);

%framewise displacement
motion_d=motion.*0;
motion_d(2:end,:)=diff(motion);
FD=sum(abs(motion_d),2);

%square of motion params
motion_square=motion.*motion;

%motion t-1
motion_pre_t=motion.*0;
motion_pre_t(2:end,:)=motion(1:end-1,:);

%motion t-1 squared
motion_pre_t_square=motion_pre_t.*motion_pre_t;

%creating censor file: bad TRs and the next TR are also bad
censor=ones(size(FD));
bad=find(FD>=proj.param.mri.FD_thresh);
censor(bad)=0;
if(numel(bad)>0)
    if(bad(end)==numel(censor)) %make sure not to increase length
                                %of censor file ***CHANGE from earlier***
        censor(bad(1:end-1)+1)=0; 
    else
        censor(bad+1)=0;
    end
end

%find single good TRs and censor them too
f=find(censor==1);
f_diff=f.*0;
f_diff(2:end)=diff(f);
bad_fs=[];
if isempty(f)==0
    for bad_loop=1:numel(f)
        if bad_loop~=numel(f)
            if f_diff(bad_loop)~=1&&f_diff(bad_loop+1)~=1
                bad_fs=[bad_fs f(bad_loop)];
            end
        elseif bad_loop==numel(f)
            if f_diff(bad_loop)~=1
                bad_fs=[bad_fs f(bad_loop)];
            end
        end
    end
end

censor(bad_fs)=0;

%% This is old way.  Better way is to use FD more wisely
%% during analysis.
%% %censor entire run if less than 50% of data is usable
%% run_length=numel(censor);
%% if numel(find(censor==1))./numel(censor)<.5
%%    censor(1:end)=0;
%% end

num_TRs=sum(censor);

motion_path=[proj.path.mri.mri_clean,study,'_',subj,'/',task, ...
            '/',scan,'/',study,'.',subj,'.',task,'.',scan];

dlmwrite([motion_path '.FD.1D'],FD,' ');
dlmwrite([motion_path '.motion_derivative.1D'],motion_d,' ');
dlmwrite([motion_path '.censor.1D'],censor,' ');
dlmwrite([motion_path '.num_TRs'],num_TRs,' ');
dlmwrite([motion_path '.motion.square.1D'],motion_square,' ');
dlmwrite([motion_path '.motion_pre_t.1D'],motion_pre_t,' ');
dlmwrite([motion_path '.motion_pre_t_square.1D'],motion_pre_t_square,' ');