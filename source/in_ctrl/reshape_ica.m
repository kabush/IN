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
%% SENSORY/MOTOR: ICs 6-12
%% HIGHER COGNITION: ICs 13-18
%%
%% We have ACC related regions in the ICAs that we want to
%% remove. (clusterizing with NN=1, voxels=2)
%% ICA2,  slight intersection (subtr. removes a few voxels): SUBTR OK
%% ICA4,  cluster #4  (1244 voxels):                         REMOVE CLUSTER
%% ICA5,  cluster #5:                                :       SUBTR OK
%% ICA6,  cluster #1  (cut out a big chunk of SMA???):       REMOVE CLUSTER
%% ICA7,  slight intersection (subtr. removes a few voxels): SUBTR OK
%% ICA13, intersection (subtr. removes dozens of voxels):    SUBTR OK
%% ICA15, cluster #'s, 12 & 16;                              REMOVE 2xCLUSTER
%% ICA17, cluster #20  (cut out a big chunk):                REMOVE CLUSTER
%% ICA18, cluster #4 (cut out a big chunk of pre-sMA???):    REMOVE CLUSTER

clust_rmv = [-1,-1,-1,4,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,12,-1,20,4];
clust_xtra = 15;
clust_rmv_xtra = 16;

ic_seq=1:18;

for iseq=1:numel(ic_seq)

    i = ic_seq(iseq);
    clust_rmv_i = clust_rmv(iseq);

    %% orient ICA to RAI like the rest of the data
    eval(['! 3dresample -orient RAI -prefix ',proj.path.code,'tmp/' ...
                        'orient_thresh_zstatd20_',num2str(i),'.nii.gz ' ...
                        '-input ',proj.path.code,'tmp/thresh_zstat',num2str(i),'.nii.gz ']);
    
    %% change from 1x1x1 to 3x3x3 voxel sizes
    eval(['! 3dfractionize -template ',proj.path.mri.gm_mask,'group_gm_mask.nii ' ...
          '-input ',proj.path.code,'tmp/orient_thresh_zstatd20_',num2str(i),'.nii.gz ' ...
          '-prefix ',proj.path.code,'tmp/' ...
          'orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz -clip .2']);

    if(clust_rmv_i<0)

        %%Just substract ACC from mask.

        %% Create Full-ICA Mask
        eval(['! 3dcalc -a ',proj.path.code,'tmp/orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz ' ...
              '-expr  ''bool(a)'' -prefix ' ...
              ,proj.path.code,'tmp/full_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz']);
        
        %% Subtract out dACC from Mask
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
              '-expr  ''equals(a,',num2str(clust_rmv_i),')'' -prefix ' ...
              ,proj.path.code,'tmp/roi_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz']);            

        %% Create ICA Mask without dACC 
        eval(['! 3dcalc -a ',proj.path.code,'tmp/full_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz ' ...
              '-b ',proj.path.code,'tmp/roi_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz ' ...
              '-expr  ''a-b'' -prefix ' ...
              ,proj.path.code,'tmp/sng_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz']);

        if(i==clust_xtra) 

            %% *** Because this IC has two small ROIs within the
            %% dACC we have to remove the 2nd cluster 

            eval(['! mv ',proj.path.code,'tmp/sng_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz ' ...
                 ,proj.path.code,'tmp/xsng_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz']);

            %% Create dACC(or partial dACC)-ICA  Mask
            eval(['! 3dcalc -a ',proj.path.code,'tmp/clst_orient_thresh_zstatd20_',num2str(i),'_3x3x3.nii.gz ' ...
                  '-expr  ''equals(a,',num2str(clust_rmv_xtra),')'' -prefix ' ...
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