function Step6A_loop_ves_ntec_without_meanint(czi_image_folder, image_name,...
    pre_processed_dir_name, results_file_name_px, results_file_name_um, vessel_processed_im,...
    org_particle_image, neighbourhood_dir_name, vessel_analysis_dir_name, ...
    metadata_file_name)

if ~exist('neighbourhood_dir_name','var')
    neighbourhood_dir_name = 'Neighbourhood_Analysis';
end
neighbourhood_dir_name

% 

if ~exist ('vessel_analysis_dir_name', 'var')
    vessel_analysis_dir_name = 'Vessel_Analysis';
end

if ~exist ('metadata_file_name', 'var')
    metadata_file_name = '_iso_info.csv';
end

image_seg_folder = czi_image_folder;
cd (image_seg_folder)
datetime()
    display (['Processing Vessels ' image_name])

    image_folder_dir = strcat(image_seg_folder,'\', image_name);
    neigh_dir = strcat(image_folder_dir, '\', neighbourhood_dir_name);
    save_dir = strcat(neigh_dir, '\', vessel_analysis_dir_name);
    
    pre_processed_dir = strcat(image_folder_dir, '\', pre_processed_dir_name);
    
    if exist(save_dir, 'dir')~=7
        mkdir(save_dir);
    end
    
    completed_results_path = strcat(save_dir, '\', image_name, results_file_name_px);
    
    if exist(completed_results_path, 'file') == 2
        display (['Vessel Analysis Already Complete ']);
    else
    


        cd(image_folder_dir)

        particles = org_particle_image;
        metadata_name = strcat(image_name,metadata_file_name);
        metadata = readtable(metadata_name,delimitedTextImportOptions);
        px_per_um = str2num(metadata.ExtraVar1{2})*1E6;

        cd(pre_processed_dir)
        vessel_image = vessel_processed_im;



       STEP6B_vessel_feature_analysis_func_without_meanint(vessel_image, particles, px_per_um, save_dir, image_name, results_file_name_px, results_file_name_um );
    end 
       datetime()
       
end