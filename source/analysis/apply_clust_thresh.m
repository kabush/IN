%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

function [] = apply_clust_thresh(proj,affect_name)

var_names = {'err','cnf','evc','pel','pro'};

for i=1:numel(var_names)

    var_name = var_names{i};
    disp(var_name);

    %eval(['! cd ',proj.path.analysis.in_clust_thresh]);
    eval(['! tail -n 1 ',proj.path.analysis.in_clust_thresh,'clust_size_', ...
                        affect_name,'_',var_name,'.NN1_bisided.1D > ./tmp/tmp.txt']);
    clust_size_array = dlmread(['./tmp/tmp.txt'],' ',0,1);
    clust_size = clust_size_array(5); %fifth value to be parsed out.

    %% Cluster threshold 3dlme findins

    cmd = ['! 3dClusterize -nosum -1Dformat -inset ', ...
           proj.path.analysis.in_3dlme, ...
           'lme_',affect_name,'_',var_name,'+tlrc -idat 4 -ithr 4 -NN 1 -clust_nvox ', ...
           num2str(clust_size),' -bisided -3.287 3.287 -pref_map ', ...
           proj.path.analysis.in_clust_thresh,'clust_mask_',affect_name, ...
           '_',var_name,' -pref_dat ', ...
           proj.path.analysis.in_clust_thresh,'clust_map_', ...
           affect_name,'_',var_name,' > ',proj.path.analysis.in_clust_thresh,'clust_info_',...
           affect_name,'_',var_name,'.1D'];

    disp(cmd)
    eval(cmd)

    % Transfer masks to *.nii and move
    cmd = ['! 3dAFNItoNIFTI ', ...
           proj.path.analysis.in_clust_thresh,'clust_mask_',affect_name, ...
           '_',var_name,'+tlrc'];
    disp(cmd);
    eval(cmd);
    
    cmd = ['! mv ./clust_mask_',affect_name,'_',var_name,...
           '.nii ',proj.path.analysis.in_clust_thresh];
    disp(cmd);
    eval(cmd);
 
    % Transfer maps to *.nii and move
    cmd = ['! 3dAFNItoNIFTI ', ...
           proj.path.analysis.in_clust_thresh,'clust_map_',affect_name, ...
           '_',var_name,'+tlrc'];
    disp(cmd);
    eval(cmd);
    
    cmd = ['! mv ./clust_map_',affect_name,'_',var_name,...
           '.nii ',proj.path.analysis.in_clust_thresh];
    disp(cmd);
    eval(cmd);

end
