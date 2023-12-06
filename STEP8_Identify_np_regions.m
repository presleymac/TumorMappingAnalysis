%Process NP channel within 100um (diameter) of nanoparticle center 
function STEP8_Identify_np_regions(image_folder_path, image_name, ...
     post_processed_dir_name, np_region_dir_name, inverted_tissue_image, ...
    mask_np_image, radius_um, pixel_size )


    
% run('C:\Users\admin\Documents\MATLAB\DIPimage 2.9\dipstart.m')

    display (['Identifying NP regions ' image_name])
    image_folder_dir = strcat(image_folder_path,'\', image_name);
    cd(image_folder_dir)
    
    post_processing_dir=(strcat(image_folder_dir,'\', post_processed_dir_name));
    np_region_dir = strcat(image_folder_dir, '\', 'Neighbourhood_Analysis_2023', '\', np_region_dir_name);
    
    if exist(np_region_dir, 'dir')~=7
        mkdir(np_region_dir)
    end 
    
    tissue_im = inverted_tissue_image;
    np_im = mask_np_image;
    
    
    np_thresh = np_im>0;
    dist_transform_np = bwdist(np_thresh);
    dist_np_rad = dist_transform_np<(radius_um/pixel_size);

    dist_transform_in_tissue = (uint16(dist_np_rad)).*tissue_im;
    dist_transform_in_tissue = dist_transform_in_tissue>0;
    
    
%     grow_np_regions = im2mat(dip_growregions(uint32(np_thresh),[],dist_transform_in_tissue,3,100,'low_first'));
    

    save_path =  strcat(np_region_dir, '\', image_name, '_', num2str(radius_um), '_um_radius_np_regions_2023.tiff');
    crop_path = strcat(np_region_dir, '\', image_name, '_',  num2str(radius_um), 'um_radius_np_regions_tissue_cropped_2023.tiff');
    
    

    full_size_name= save_path;
           clear options;
           options.overwrite = true;
           options.compress = 'lzw';
           saveastiff(uint16(dist_np_rad), full_size_name, options);              

     cropped_name= crop_path;
           clear options;
           options.overwrite = true;
           options.compress = 'lzw';
           saveastiff(uint16(dist_transform_in_tissue), cropped_name, options);
% 
%     
    
    
%     break
