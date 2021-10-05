function [dcmp,indx] = decompose_in(proj,label_id,hd)

%% ----------------------------------------
%% Calculate dynamcs from IN

%% Specialized IN ids
in_ids = find(label_id==1);
cue_ids = find(label_id==2);
feel_ids = find(label_id==3);
rest_ids = find(label_id==4);

%% Length params
% Nstim = 1;
% Ncue = 1;
% Nfeel = 4;
% Nrest = 1;
Nstim = proj.param.trg.cogdyn.n_stim;
Ncue = proj.param.trg.cogdyn.n_cue;
Nfeel = proj.param.trg.cogdyn.n_feel;
Nrest = proj.param.trg.cogdyn.n_rest;
Ntot = Nstim+Ncue+Nfeel+Nrest;

dcmp = struct();    
indx = struct();

%%----------------------------------------
%% Construct IN individual pieces 

%%Get stim
dcmp.stim = hd(in_ids);
indx.stim = in_ids;

%%Get cue
dcmp.cue = hd(cue_ids);
indx.cue = cue_ids;

%%Get ctrl response
tmp_feel = hd(feel_ids);
dcmp.feel = reshape(tmp_feel,Nfeel,numel(feel_ids)/Nfeel)';
indx.feel = reshape(feel_ids,Nfeel,numel(feel_ids)/Nfeel)';

%%Get rest
dcmp.rest = hd(rest_ids);
indx.rest = rest_ids;

%%----------------------------------------
%% Construct IN trajectories
dcmp.h = [];
indx.h = [];
for i=1:numel(dcmp.stim)
    traj = [dcmp.stim(i),dcmp.cue(i),dcmp.feel(i,:),dcmp.rest(i)];
    traj_idx = [indx.stim(i),indx.cue(i),indx.feel(i,:),indx.rest(i)];
    dcmp.h = [dcmp.h;traj];
    indx.h = [indx.h;traj_idx];
end

%%Construct plant derivative
dcmp.dh = zeros(numel(in_ids),numel(2:(Ntot-1)));
indx.dh = zeros(numel(in_ids),numel(2:(Ntot-1)));
for i = 1:numel(in_ids)
    for j = 2:(Ntot-1)
        dcmp.dh(i,j-1) = (dcmp.h(i,j+1)-dcmp.h(i,j-1))/2;
        indx.dh(i,j-1) = indx.h(i,j);
    end
end

%%Construct plant 2nd derivative
dcmp.d2h = zeros(numel(in_ids),numel(3:(Ntot-2)));
indx.d2h = zeros(numel(in_ids),numel(3:(Ntot-2)));
for i = 1:numel(in_ids)
    for j = 3:(Ntot-2)
        dcmp.d2h(i,j-2) = (dcmp.dh(i,j)-dcmp.dh(i,j-2))/2;
        indx.d2h(i,j-2) = indx.dh(i,j-1);
    end
end

%%Construct plant 3rd derivative
dcmp.d3h = zeros(size(dcmp.h,1),numel(4:(Ntot-3)));
indx.d3h = zeros(size(dcmp.h,1),numel(4:(Ntot-3)));
for i = 1:numel(in_ids)
    for j = 4:(Ntot-3)
        dcmp.d3h(i,j-3) = (dcmp.d2h(i,j-1)-dcmp.d2h(i,j-3))/2;
        indx.d3h(i,j-3) = indx.d2h(i,j-2);
    end
end

%%----------------------------------------
%%Construct IN-based error trajectories
dcmp.err = 0*dcmp.h;
indx.err = indx.h;
for i=1:numel(in_ids)
    dcmp.err(i,:) = dcmp.stim(i)-dcmp.h(i,:);
end

%%Construct error derivative
dcmp.derr = 0*dcmp.dh;
indx.derr = 0*dcmp.dh;
for i=1:numel(in_ids)
    for j = 2:(Ntot-1)
        dcmp.derr(i,j-1) = (dcmp.err(i,j+1)-dcmp.err(i,j-1))/2;
        indx.derr(i,j-1) = indx.err(i,j);
    end
end

%%Construct error 2nd derivative
dcmp.d2err = 0*dcmp.d2h;
indx.d2err = 0*dcmp.d2h;
for i=1:numel(in_ids)
    for j = 3:(Ntot-2)
        dcmp.d2err(i,j-2) = (dcmp.derr(i,j)-dcmp.derr(i,j-2))/2;
        indx.d2err(i,j-1) = indx.derr(i,j);
    end
end

%%Construct error 3rd derivative
dcmp.d3err = 0*dcmp.d3h;
indx.d3err = 0*dcmp.d3h;
for i=1:numel(in_ids)
    for j = 4:(Ntot-3)
        dcmp.d3err(i,j-3) = (dcmp.d2err(i,j-1)-dcmp.d2err(i,j-3))/2;
        indx.d3err(i,j-1) = indx.d2err(i,j);
    end
end
