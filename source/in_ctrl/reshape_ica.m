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

logger(['*******************************************'],proj.path.logfile);
logger([' Reshape Ray 2013 ACC ICA to match our data'],proj.path.logfile);
logger(['*******************************************'],proj.path.logfile);

%% ----------------------------------------
%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.ctrl.in_ica]);
    eval(['! rm -rf ',proj.path.ctrl.in_ica]);
    disp(['Creating ',proj.path.ctrl.in_ica]);
    eval(['! mkdir ',proj.path.ctrl.in_ica]);
end

%% ----------------------------------------
%% ----------------------------------------
%% Configure the ACC mask

%% copy ica to tmp
eval(['! cp ',proj.path.atlas,'ray2013/ica70/maps/' ...
                    'thresh_zstatd70_17.nii.gz ',proj.path.code,'tmp/']);
eval(['! cp ',proj.path.atlas,'TT/TT_icbm452_orig.nii ',proj.path.code,'tmp/']);

%% rotate the ICA to match the rotation of our fMRI data
eval(['! 3dresample -orient RAI -prefix ',proj.path.code,'tmp/orient_thresh_zstatd70_17.nii.gz ' ...
      '-input ',proj.path.code,'tmp/thresh_zstatd70_17.nii.gz ']);

%% find ica values greater than threshold to achieve largest single
%% ROI (>5.6734)
eval(['! 3dcalc -a ',proj.path.code,'tmp/orient_thresh_zstatd70_17.nii.gz ' ...
                    '-expr  ''ispositive(a-5.6734)'' -prefix ' ...
                    ,proj.path.code,'tmp/frc_orient_thresh_zstatd70_17.nii.gz']);

%% change from 1x1x1 to 3x3x3 voxel sizes
eval(['! 3dfractionize -template ',proj.path.mri.gm_mask,'group_gm_mask.nii ' ...
      '-input ',proj.path.code,'tmp/frc_orient_thresh_zstatd70_17.nii.gz ' ...
      '-prefix ',proj.path.code,'tmp/' ...
      'frc_orient_thresh_zstatd70_17_3x3x3.nii.gz -clip .2']);

%% create mask 
eval(['! 3dcalc -a ',proj.path.code,'tmp/frc_orient_thresh_zstatd70_17_3x3x3.nii.gz ' ...
      '-expr  ''bool(a)'' -prefix ' ...
      ,proj.path.code,'tmp/sng_orient_thresh_zstatd70_17_3x3x3.nii.gz']);

%% clusterize to remove everything but the dACC
eval(['! 3dClusterize -nosum -1Dformat -inset ',proj.path.code, ...
      'tmp/sng_orient_thresh_zstatd70_17_3x3x3.nii.gz ' ...
      '-idat 0 -ithr 0 -NN 1 -clust_nvox 100 -1sided RIGHT_TAIL 0.5 ' ...
      '-pref_map ',proj.path.code,'tmp/clst_sng_orient_thresh_zstatd70_17_3x3x3.nii.gz']);


%% ----------------------------------------
%% ----------------------------------------
%% Configure the State mask (for EVC)

%% copy ica to tmp
eval(['! cp ',proj.path.atlas,'ray2013/ica20/maps/*.gz ',proj.path.code,'tmp/']);

%% rotate the ICA to match the rotation of our fMRI data

%% ----------------------------------------
%% As per Ray et al., (2013)
%% EMOTION/INTEROCEPTION REGIONS: ICs 1-5
%% HIGHER COGNITION: ICs 13-18
%%
%% We have ACC related regions in the ICAs that we want to
%% remove. (clusterizing with NN=1, voxels=1)
%% ICA2, cluster #14 (46 voxels)
%% ICA4, cluster #4  (1244 voxels)
%% ICA5, cluster #6  (137 voxels)
ic_seq=[1:5,13:18];
clust_seq = [-1,14,-1,4,6,-1,-1,12,-1,20,4];  %%IC=15 also removes 16
xclust_id = 16;

for iseq=1:numel(ic_seq)

    i = ic_seq(iseq);
    clust_i = clust_seq(iseq);

    %% orient ICA to RAI like the rest of the data
    eval(['! 3dresample -orient RAI -prefix ',proj.path.code,'tmp/' ...
                        'orient_thresh_zstatd20_',num2str(i),'.nii.gz ' ...
                        '-input ',proj.path.code,'tmp/thresh_zstat',num2str(i),'.nii.gz ']);
    
    %% change from 1x1x1 to 3x3x3 voxel sizes
    eval(['! 3dfractionize -template ',proj.path.mri.gm_mask,'group_gm_mask.nii ' ...
          '-input ',proj.path.code,'tmp/orient_thresh_zstatd20_',num2str(i),'.nii.gz ' ...
          '-prefix ',proj.path.code,'tmp/' ...
          'orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz -clip .2']);

    if(i==13)

        %% **** Because the overlap is just a few voxels out of
        %% many hundreds, we subtract out the dACC voxels from this
        %% IC rather than remove the ROI of the IC that is
        %% intersecting dACC.  ***

        %% Create Full-ICA Mask
        eval(['! 3dcalc -a ',proj.path.code,'tmp/orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz ' ...
              '-expr  ''bool(a)'' -prefix ' ...
              ,proj.path.code,'tmp/full_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz']);
  
        %% Subject out dACC from Mask
        eval(['! 3dcalc -a ',proj.path.code,'tmp/full_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz ' ...
              '-b ',proj.path.code,'tmp/clst_sng_orient_thresh_zstatd70_17_3x3x3.nii.gz ' ...
              '-expr  ''a-b'' -prefix ' ...
              ,proj.path.code,'tmp/diff_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz']);

        %% Keep only the positive
        eval(['! 3dcalc -a ',proj.path.code,'tmp/diff_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz ' ...
              '-expr  ''ispositive(a)'' -prefix ' ...
              ,proj.path.code,'tmp/sng_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz']);

    else

        %% clusterize to remove everything but the dACC
        eval(['! 3dClusterize -nosum -1Dformat -inset ',proj.path.code, ...
              'tmp/orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz ' ...
              '-idat 0 -ithr 0 -NN 1 -clust_nvox 2 -1sided RIGHT_TAIL 0.5 ' ...
              '-pref_map ',proj.path.code,'tmp/clst_orient_thresh_zstatd20_', ...
              num2str(i),'_3x3x3.nii.gz']);
        
        %% Create Full-ICA Mask
        eval(['! 3dcalc -a ',proj.path.code,'tmp/clst_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz ' ...
              '-expr  ''bool(a)'' -prefix ' ...
              ,proj.path.code,'tmp/full_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz']);
        
        %% Create dACC(or partial dACC)-ICA  Mask
        eval(['! 3dcalc -a ',proj.path.code,'tmp/clst_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz ' ...
              '-expr  ''equals(a,',num2str(clust_i),')'' -prefix ' ...
              ,proj.path.code,'tmp/roi_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz']);            

        %% Create ICA Mask without dACC 
        eval(['! 3dcalc -a ',proj.path.code,'tmp/full_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz ' ...
              '-b ',proj.path.code,'tmp/roi_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz ' ...
              '-expr  ''a-b'' -prefix ' ...
              ,proj.path.code,'tmp/sng_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz']);

        if(i==15) 

            %% *** Because this IC has two small ROIs within the
            %% dACC we have to remove the 2nd ***6

            eval(['! mv ',proj.path.code,'tmp/sng_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz ' ...
                 ,proj.path.code,'tmp/xsng_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz']);

            %% Create dACC(or partial dACC)-ICA  Mask
            eval(['! 3dcalc -a ',proj.path.code,'tmp/clst_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz ' ...
                  '-expr  ''equals(a,',num2str(xclust_id),')'' -prefix ' ...
                  ,proj.path.code,'tmp/xroi_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz']);            

            %% Create ICA Mask without dACC 
            eval(['! 3dcalc -a ',proj.path.code,'tmp/xsng_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz ' ...
                  '-b ',proj.path.code,'tmp/xroi_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz ' ...
                  '-expr  ''a-b'' -prefix ' ...
                  ,proj.path.code,'tmp/sng_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz']);

        end

    end
    
    %% move to permanent storage
    eval(['! mv ',proj.path.code,'tmp/sng_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz ', ... 
          proj.path.ctrl.in_ica]);
    
end

%% move to permanent storage
eval(['! mv ',proj.path.code,'tmp/clst_sng_orient_thresh_zstatd70_17_3x3x3.nii.gz ',proj.path.ctrl.in_ica]);
eval(['! mv ',proj.path.code,'tmp/TT_icbm452_orig.nii ',proj.path.ctrl.in_ica]);

%% clean-up
eval(['! rm ',proj.path.code,'tmp/*']);