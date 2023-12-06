function [img_skel, branch_pts_image, labeled_vess_segs, vess_radius, ...
    binary_hs_image, hotspots_pts_image_dil, results_table_um_new_table,...
    hs_diameter_table] = STEP6B_vessel_feature_analysis_func_without_meanint(processed_ves,...
    original_particle_image, px_per_um, save_dir,sample_name, results_name_px, results_name_um)


%%% Code orginally designed by Kingston et al., 2021. Origial code can 
%%% be found at:
%%% https://github.com/BenKingston/NTEC_Vessel_Analysis/



tic
shortfile = sample_name;
display(['Analyzing vessels for ' shortfile])


%%%% Vessel skeletonization and vessel chopping into
%%%% individual vessel segments

ves_image = processed_ves;
ves_thresh = thresh(ves_image, 1);

slow_img_skel = bwskel(ves_thresh,'MinBranchLength',15);
slow_img_skel = slow_img_skel>0;

img_skel = slow_img_skel;

branch_pts = bwmorph3(slow_img_skel,'branchpoints'); 
vessel_seg_skel = slow_img_skel-branch_pts;
end_pts = bwmorph3(vessel_seg_skel,'endpoints');
label_vessel_seg_skel = bwlabeln(vessel_seg_skel);

stats = regionprops3(label_vessel_seg_skel,'Volume','VoxelIdxList');

labeled_vess_segs = im2mat(dip_growregions(uint32(label_vessel_seg_skel),[],ves_thresh,3,10,'low_first'));

labeled_vess_seg_ext_regions = im2mat(dip_growregions(uint32(labeled_vess_segs),[],[],3,2,'low_first'));

SE = strel('sphere',6);           
branch_pts_image = imdilate(branch_pts, SE);

%%Hotspot assignment - could be removed as not used in further analysis
%%kept to be consitant with the original analysis
particles = original_particle_image;

dist_tform_vessel = bwdist(ves_thresh);

dist_threshold = dist_tform_vessel<5;
np_near_vess = particles.*uint16(dist_threshold);

vess_NP_Pixels = np_near_vess(np_near_vess > 0);
threshold_NP_near_ves = graythresh(vess_NP_Pixels);
treshold_np = imbinarize(np_near_vess,threshold_NP_near_ves*3);
label_hotspots = bwlabeln(treshold_np);

hotspots_analysis = regionprops3(label_hotspots,np_near_vess,'WeightedCentroid','EquivDiameter');
binary_hs_image = zeros(size(np_near_vess));

rounded_pts = round(hotspots_analysis.WeightedCentroid);

for r = 1:size(rounded_pts,1)
binary_hs_image(rounded_pts(r,2),rounded_pts(r,1),rounded_pts(r,3))=1;
end

label_hs = bwlabeln(binary_hs_image);
SE_hs = strel('sphere',5);           
hotspots_pts_image_dil = imdilate(label_hs, SE_hs);

hotspot_ves_segs = regionprops3(labeled_vess_seg_ext_regions,binary_hs_image,'VoxelValues');
array_hotspot_ves_segs = table2array(hotspot_ves_segs);

hs_per_ves_seg = zeros(size(array_hotspot_ves_segs,1),1);

for t=1:length(array_hotspot_ves_segs)
hs_per_ves_seg(t)=sum(array_hotspot_ves_segs{t},1);
end

hs_per_ves_seg_table = array2table(hs_per_ves_seg);

hs_per_ves_seg_table.Properties.VariableNames = {'Num_hotsposts_per_ves_seg'};

%%Assign pixel dimensions in um
px_per_um = px_per_um;

%%hotspot diameter
HS_diameter_table = hotspots_analysis.EquivDiameter;

results_table_dia_um_new = [HS_diameter_table(:,1).*px_per_um];
results_table_dia_um_new = array2table(results_table_dia_um_new);

results_table_dia_um_new.Properties.VariableNames = {'Hotspot_diameter_um'};
hs_diameter_table = results_table_dia_um_new;

%%Vess diameter distance transform
vess_neg = ves_thresh~=1;
vess_neg_dia = bwdist(vess_neg);

vess_radius = vess_neg_dia; 

%%%%Vessel Length using skeleton where length = volume
ves_length = regionprops3(label_vessel_seg_skel,'Volume');

ves_length.Properties.VariableNames = {'Vessel_length_px'};

%%%%Vessel Diameter from average radius mulitplied by 2
ves_rad = regionprops3(label_vessel_seg_skel,vess_neg_dia,'MeanIntensity');
ves_dia_array = (table2array(ves_rad))*2;
ves_dia_table = array2table(ves_dia_array);

ves_dia_table.Properties.VariableNames = {'Vessel_diameter_px'};

%%%%Vessel Surface area and Volume of each vessel segment
ves_sa_vol = regionprops3(labeled_vess_segs,'Volume','SurfaceArea','Centroid');
vess_sa_vol_comb = [ves_sa_vol.Volume ves_sa_vol.SurfaceArea];
vess_sa_vol_table = array2table(vess_sa_vol_comb);

vess_sa_vol_table.Properties.VariableNames = {'Vessel_Volume_px' 'Vessel_Surface_Area_px'};

%%%%Vessel Surface area:Volume
ves_sa_to_vol_ratio = ves_sa_vol.SurfaceArea./ves_sa_vol.Volume;
ves_sa_to_vol_ratio_table = array2table(ves_sa_to_vol_ratio);

ves_sa_to_vol_ratio_table.Properties.VariableNames = {'Vessel_SA_to_Vol_Ratio'};

%%%%Vessel Distance to nearest junction point 
branch_pts_Centroid = regionprops3(branch_pts,'Centroid');
dist_allpopints = pdist2(ves_sa_vol.Centroid,branch_pts_Centroid.Centroid);
ves_dist_to_branch_pt = min(dist_allpopints,[],2);
ves_dist_to_branch_pt_table = array2table(ves_dist_to_branch_pt);

ves_dist_to_branch_pt_table.Properties.VariableNames = {'Vessel_dist_to_branch_pt_px'};

%%%%Vessel Distance to nearest vessel segment, nearest 5 vessels, nearest
%%%%10 vessels 
dist_all_ves = pdist(ves_sa_vol.Centroid);
dist_all_ves_sq = squareform(dist_all_ves);

dist_all_ves_min = mink(dist_all_ves_sq,2,2);
dist_all_ves_min_5 = mink(dist_all_ves_sq,6,2);
dist_all_ves_min_10 = mink(dist_all_ves_sq,11,2);

dist_all_ves_min_5 = dist_all_ves_min_5(:,2:6);
dist_all_ves_min_10 = dist_all_ves_min_10(:,2:11);
dist_all_ves_min = dist_all_ves_min(:,2);

dist_all_ves_min_5_mean = mean(dist_all_ves_min_5,2);
dist_all_ves_min_10_mean = mean(dist_all_ves_min_10,2);

dist_all_ves_closest_table = array2table(dist_all_ves_min);
dist_all_ves_min_5_mean_table = array2table(dist_all_ves_min_5_mean);
dist_all_ves_min_10_mean_table = array2table(dist_all_ves_min_10_mean);

dist_all_ves_closest_table.Properties.VariableNames = {'Vessel_dist_to_nearest_ves_px'};
dist_all_ves_min_5_mean_table.Properties.VariableNames = {'Vessel_dist_to_nearest_5_ves_px'};
dist_all_ves_min_10_mean_table.Properties.VariableNames = {'Vessel_dist_to_nearest_10_ves_px'};

%%%%Create final table

results_table = [ves_length ves_dia_table vess_sa_vol_table ves_sa_to_vol_ratio_table ves_dist_to_branch_pt_table dist_all_ves_closest_table dist_all_ves_min_5_mean_table dist_all_ves_min_10_mean_table hs_per_ves_seg_table];

results_table_um = table2array(results_table);

results_table_um_new = [results_table_um(:,1).*px_per_um results_table_um(:,2).*px_per_um results_table_um(:,3).*(px_per_um*px_per_um*px_per_um) results_table_um(:,4).*(px_per_um*px_per_um) results_table_um(:,5) results_table_um(:,6).*px_per_um results_table_um(:,7).*px_per_um results_table_um(:,8).*px_per_um  results_table_um(:,9) results_table_um(:,10) ];

results_table_um_new_table = array2table(results_table_um_new);


results_table_um_new_table.Properties.VariableNames = {'Vessel_length_um' 'Vessel_diameter_um' 'Vessel_vol_um' 'Vessel_SA_um' 'Vessel_SA_to_Vol_ratio' 'Vessel_dist_to_branch_pt_um' 'Vessel_dist_to_nearest_vessel_um' 'Vessel_dist_to_nearest_5_vessel_um' 'Vessel_dist_to_nearest_10_vessel_um' 'Num_hotspots_per_ves_seg'};


%%Save Tables in save_dir

cd(save_dir)
table_name = strcat(shortfile,results_name_px);
table_name_um = strcat(shortfile,results_name_um);

writetable(results_table,table_name);
writetable(results_table_um_new_table,table_name_um);

table_name_dia = strcat(shortfile,'_Results_hotspot_dia_analysis.xlsx');

writetable(hs_diameter_table,table_name_dia);

% table_name_diff = strcat(shortfile,'_Diffusion_data.csv');
% writetable(result_test_table_diff,table_name_diff);

%%%% Write Image files to save_dir

cd(save_dir)

img_skel_name = strcat(shortfile,'_skeleton.tif');
clear options;
            options.overwrite = true;
            options.compress = 'lzw';
            saveastiff(uint16(slow_img_skel), img_skel_name, options);
            
branch_pts_image_name = strcat(shortfile,'_branch_pts.tif');
clear options;
            options.overwrite = true;
            options.compress = 'lzw';
            saveastiff(uint16(branch_pts_image), branch_pts_image_name, options);
            
labeled_vess_segs_name = strcat(shortfile,'_labeled_vess_segs.tif');
clear options;
            options.overwrite = true;
            options.compress = 'lzw';
            saveastiff(uint16(labeled_vess_segs), labeled_vess_segs_name, options);            

vess_radius_segs_name = strcat(shortfile,'_vess_radius.tif');
clear options;
            options.overwrite = true;
            options.compress = 'lzw';
            saveastiff(uint16(vess_neg_dia), vess_radius_segs_name, options);   

vess_hs_segs_name = strcat(shortfile,'_bin_hs_img.tif');
clear options;
            options.overwrite = true;
            options.compress = 'lzw';
            saveastiff(uint16(binary_hs_image), vess_hs_segs_name, options);  

vess_dil_hs_segs_name = strcat(shortfile,'_dilated_hs_img.tif');
clear options;
            options.overwrite = true;
            options.compress = 'lzw';
            saveastiff(uint16(hotspots_pts_image_dil), vess_dil_hs_segs_name, options);   
            
                  
toc
end