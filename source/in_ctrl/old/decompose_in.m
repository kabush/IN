function [dcmp] = decompose_in(proj,label_id,hd)

%% ----------------------------------------
%% Calculate dynamcs from IN

%% Specialized IN idss
in_ids = find(label_id==1);
cue_ids = find(label_id==2);
feel_ids = find(label_id==3);
rest_ids = find(label_id==4);

%% Length params
Nstim = 1;
Ncue = 1;
Nfeel = 4;
Nrest = 1;
Ntot = Nstim+Ncue+Nfeel+Nrest;

dcmp = struct();    

%%----------------------------------------
%% Construct IN individual pieces 

%%Get stim
dcmp.stim = hd(in_ids);

%%Get cue
dcmp.cue = hd(cue_ids);

%%Get ctrl response
tmp_feel = hd(feel_ids);
dcmp.feel = reshape(tmp_feel,Nfeel,numel(feel_ids)/Nfeel)';

%%Get rest
dcmp.rest = hd(rest_ids);

%%----------------------------------------
%% Construct IN trajectories
dcmp.h = [];
for i=1:numel(dcmp.stim)
    traj = [dcmp.stim(i),dcmp.cue(i),dcmp.feel(i,:),dcmp.rest(i)];
    dcmp.h = [dcmp.h;traj];
end

%%Construct plant derivative
dcmp.dh = zeros(numel(in_ids),numel(2:(Ntot-1)));
for i = 1:numel(in_ids)
    for j = 2:(Ntot-1)
        dcmp.dh(i,j-1) = (dcmp.h(i,j+1)-dcmp.h(i,j-1))/2;
    end
end

%%Construct plant 2nd derivative
dcmp.d2h = zeros(numel(in_ids),numel(3:(Ntot-2)));
for i = 1:numel(in_ids)
    for j = 3:(Ntot-2)
        dcmp.d2h(i,j-2) = (dcmp.dh(i,j)-dcmp.dh(i,j-2))/2;
    end
end

%%Construct plant 3rd derivative
dcmp.d3h = zeros(size(dcmp.h,1),numel(4:(Ntot-3)));
for i = 1:numel(in_ids)
    for j = 4:(Ntot-3)
        dcmp.d3h(i,j-3) = (dcmp.d2h(i,j-1)-dcmp.d2h(i,j-3))/2;
    end
end

%%----------------------------------------
%%Construct IN-based error trajectories
dcmp.err = 0*dcmp.h;
for i=1:numel(in_ids)
    dcmp.err(i,:) = dcmp.stim(i)-dcmp.h(i,:);
end

%%Construct error derivative
dcmp.derr = 0*dcmp.dh;
for i=1:numel(in_ids)
    for j = 2:(Ntot-1)
        dcmp.derr(i,j-1) = (dcmp.err(i,j+1)-dcmp.err(i,j-1))/2;
    end
end

%%Construct error 2nd derivative
dcmp.d2err = 0*dcmp.d2h;
for i=1:numel(in_ids)
    for j = 3:(Ntot-2)
        dcmp.d2err(i,j-2) = (dcmp.derr(i,j)-dcmp.derr(i,j-2))/2;
    end
end

%%Construct error 3rd derivative
dcmp.d3err = 0*dcmp.d3h;
for i=1:numel(in_ids)
    for j = 4:(Ntot-3)
        dcmp.d3err(i,j-3) = (dcmp.d2err(i,j-1)-dcmp.d2err(i,j-3))/2;
    end
end
