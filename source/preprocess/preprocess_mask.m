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

%% Set-up Directory Structure
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.mri.gm_mask]);
    eval(['! rm -rf ',proj.path.mri.gm_mask]);
    disp(['Creating ',proj.path.mri.gm_mask]);
    eval(['! mkdir ',proj.path.mri.gm_mask]);
end

%% Create the subjects to be analyzed (possible multiple studies)
subjs = load_subjs(proj);

logger(['****************************************'],proj.path.logfile);
logger(['Processing GM Mask of ',num2str(numel(subjs)),' subjects'],proj.path.logfile);
logger(['****************************************'],proj.path.logfile);

%%========================================
%% Convert AFNI to NIFTI
%%========================================
for i=1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    %% debug
    logger([subj_study,':',name],proj.path.logfile);

    in_path = [proj.path.mri.mri_clean,subj_study,'_',name,'/anat/', ...
               subj_study,'.',name,'.anat.seg.fsl.MNI.GM+tlrc'];
    out_path = [proj.path.mri.gm_mask,subj_study,'.',name,'.gm.nii'];

    %% convert the afni and store in gm_mask directory
    eval(['! 3dAFNItoNIFTI ',in_path]);
    eval(['! mv ',proj.path.code,subj_study,'.',name,'.anat.seg.fsl.MNI.GM.nii ',out_path]);

end

%%========================================
%% Load Base Grey-Matter Mask (for sizing)
%%========================================
subj_study = subjs{1}.study;
name = subjs{1}.name;
base_gm_mask = load_untouch_nii([proj.path.mri.gm_mask,subj_study,'.',name,'.gm.nii']);
base_gm_mask.img = 0*base_gm_mask.img;

%%========================================
%% Sum Grey-Matter Masks for all subjects
%%========================================
for i=1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    
    %% summarize individual masks
    gm_mask = load_untouch_nii([proj.path.mri.gm_mask,subj_study,'.',name,'.gm.nii']);
    base_gm_mask.img = base_gm_mask.img + gm_mask.img;

end

%%========================================
%% Threshold by 50% subject presencce
%%========================================
Nmax = max(max(max(base_gm_mask.img)));
threshold = round(Nmax/2);
brain_size = size(base_gm_mask.img);
ids = find(base_gm_mask.img>=threshold);
final_gm_mask = base_gm_mask;
final_gm_mask.img = 0*base_gm_mask.img;
final_gm_mask.img(ind2sub(brain_size,ids)) = 1;

% %% Add important header information (important)
% Nmaps = 1;
% final_gm_mask.hdr.dime.dim=[4 brain_size Nmaps 1 1 1];
% final_gm_mask.hdr.dime.datatype=64;
% final_gm_mask.hdr.dime.bitpix=64;
% final_gm_mask.original.hdr.dime.dim=[4 brain_size Nmaps 1 1 1];
% final_gm_mask.original.hdr.dime.datatype=64;
% final_gm_mask.original.hdr.dime.bitpix=64;
% 
% %%put in RAI orientation (important)
% orient=[4 5 3];
% final_gm_mask=rri_orient(final_gm_mask,orient);%silently swear at jimmy shen

%% Test of number of GM voxels  (30K-50K)
gm_vec = vec_img_2d_nii(final_gm_mask);
logger(['Num. GM voxels: ',num2str(sum(gm_vec)),', should be 30K-50K.'],proj.path.logfile);

%% Save out grey-matter
save_untouch_nii(final_gm_mask,[proj.path.mri.gm_mask,'/group_gm_mask.nii']);

%% Indicate compeltion of this process
proj.process.mask = 1;

%% Write out amended project params
save('proj.mat');
