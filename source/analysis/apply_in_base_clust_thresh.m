%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================


function [] = apply_in_base_clust_thresh(proj,affect_name)

var_names = proj.param.ctrl.base_z_names; 
var_ids = proj.param.ctrl.base_z_ids;

for i=1:numel(var_names)

    name = var_names{i};
    id = var_ids{i};
    
    eval(['! tail -n 1 ',proj.path.analysis.in_base_clust_thresh,'clust_size_', ...
          affect_name,'_base.NN1_bisided.1D > ./tmp/tmp.txt']);
    clust_size_array = dlmread(['./tmp/tmp.txt'],' ',0,1);
    clust_size = clust_size_array(5); %fifth value to be parsed out.
    
    % Cluster threshold 3dlme findings
    cmd = ['! 3dClusterize -nosum -1Dformat -inset ', ...
           proj.path.analysis.in_base_3dlme,'lme_',affect_name, ...
           '_base+tlrc -idat ',num2str(id),' -ithr ',num2str(id), ...
           ' -NN 1 -clust_nvox ', num2str(clust_size),...
           ' -bisided -3.287 3.287 -pref_map ', ...
           proj.path.analysis.in_base_clust_thresh,'clust_mask_', ...
           affect_name, '_',name,' -pref_dat ', ...
           proj.path.analysis.in_base_clust_thresh,'clust_map_', ...
           affect_name,'_',name,' > ', ...
           proj.path.analysis.in_base_clust_thresh,'clust_info_', ...
           affect_name,'_',name,'.1D'];
    disp(cmd)
    eval(cmd)
    
    % Transfer masks to *.nii and move
    cmd = ['! 3dAFNItoNIFTI ', ...
           proj.path.analysis.in_base_clust_thresh,'clust_mask_',affect_name, ...
           '_',name,'+tlrc'];
    disp(cmd);
    eval(cmd);
    
    cmd = ['! mv ./clust_mask_',affect_name,'_',name,'.nii ',...
           proj.path.analysis.in_base_clust_thresh];
    disp(cmd);
    eval(cmd);
    
    % Transfer maps to *.nii and move
    cmd = ['! 3dAFNItoNIFTI ', ...
           proj.path.analysis.in_base_clust_thresh,'clust_map_',affect_name, ...
           '_',name,'+tlrc'];
    disp(cmd);
    eval(cmd);
    
    cmd = ['! mv ./clust_map_',affect_name,'_',name,'.nii ',...
           proj.path.analysis.in_base_clust_thresh];
    disp(cmd);
    eval(cmd);
    
end
