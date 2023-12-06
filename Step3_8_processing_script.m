%MATLAB processing script for tumour neighbourhoods runs image by image to
%completion of processing after ilastik segemntation:

%Step 3: Macrophage post processing
%Step 4: Label macrophages
% Step5: Vessel Processing
% Step6: Identify and measure vessel features
% 	6A: Identify features and extract metadata
%     6B: Run analysis script and save vessel data
% Step7: Process Nanoparticles
% Step8: Identify Nanoparticle Regions
% Step9: Identify Nanoparticle Region Centers and save images

% After this point further analysis and model building is conducted with
% Python3.6.12
run('C:\Users\admin\Documents\MATLAB\DIPimage 2.9\dipstart.m')

tic
czi_image_folder = uigetdir();
pre_processed_folder_name = 'pre_processing_2022';
post_processed_folder_name = 'post_processed_2023';
np_region_dir_name = 'NP_Region_Images';

macrophage_segmentation_end = '_iso_ch4_Simple Segmentation_2022.tiff';
macrophage_segmentation_end_downsampled = '_iso_ch4_downsampled_restored_Simple Segmentation_2022.tiff';
vessel_segmentation_end = '_iso_ch3_Simple Segmentation_2022.tiff';
org_particle_image_end = '_iso_ch5.tif';

tissue_outline_save_end = '_tissue_outline_2023.tiff';
processed_macro_str = '_post_processed_macrophages_2023.tiff';
processed_ves_str = '_post_processed_vessels_2023.tiff';
nanoparticle_image_save_str = '_post_processed_nanoparticles_2023.tiff';
inverted_tissue_outline_name = '_tissue_outline_inverted_2023.tiff';
thresholded_np_name  = '_nanoparticle_image_2x_otsu_thresholded_2023.tiff';

results_file_name_px = '_Vessel_analysis_results_pixels_2023.xlsx';
results_file_name_um = 'Vessel_analysis_results_um_2023.xlsx';

neighbourhood_dir_name = 'Neighbourhood_Analysis_2023';
vessel_analysis_dir_name = 'Vessel_Analysis_2023';
metadata_file_name = '_iso_info.csv';

radius_list = [25, 50, 100];



cd (czi_image_folder)
MSC = dir('MSC*');
UT = dir('UT*');

mergestructs = @(x,y) cell2struct([struct2cell(x),struct2cell(y)],fieldnames(x),1);
image_names = mergestructs(MSC,UT);

for image = 1:size(image_names,1)
%     tic
    [~,single_image_name] = fileparts(image_names(image).name);
    
    image_folder_dir = strcat(czi_image_folder,'\', single_image_name);
    
    cd(image_folder_dir)
    org_part_image = imreadfast(strcat(single_image_name, org_particle_image_end));
    tic
    %get paths and load images
    pre_processed_dir =  strcat(image_folder_dir, '\', pre_processed_folder_name);
    post_processed_dir =  strcat(image_folder_dir, '\', post_processed_folder_name);
    macro_seg_name = strcat(single_image_name, macrophage_segmentation_end);
    macro_seg_name_downsampled = strcat(single_image_name, macrophage_segmentation_end_downsampled);
    processed_macro_name = strcat(single_image_name, processed_macro_str);
    vessel_seg_name = strcat(single_image_name, vessel_segmentation_end);
    tissue_outline_name = strcat(single_image_name, tissue_outline_save_end);
    
    if exist(pre_processed_dir, 'dir')~=7
        mkdir(pre_processed_dir)
    end
    if exist(post_processed_dir, 'dir')~=7
        mkdir(post_processed_dir)
    end
    
    cd(pre_processed_dir)
    if exist(macro_seg_name, 'file')==2
            macro_seg_im =  imreadfast(macro_seg_name);
    elseif exist(macro_seg_name_downsampled, 'file')==2
            macro_seg_im = imreadfast(macro_seg_name_downsampled);
    else
        fprintf('No Macrophage Segmentation File')
    end
    vessel_seg_im = imreadfast(vessel_seg_name);
    
    
        metadata_name = strcat(czi_image_folder, '\', single_image_name,'\', single_image_name,metadata_file_name);
        metadata = readtable(metadata_name,delimitedTextImportOptions);
        px_per_um = str2num(metadata.ExtraVar1{2})*1E6;
    
    
    fprintf('images_loaded\n')
%     toc
%     tic
    cd(post_processed_dir)
    if exist(tissue_outline_name, 'file') ==2
        
        tissue_outline = imreadfast(tissue_outline_name);
        [macro_processed_image, tissue_outline_im]= Step3_macrophage_post_processing_function(czi_image_folder, ...
            single_image_name, pre_processed_folder_name, post_processed_folder_name, ...
            vessel_seg_im, macro_seg_im, processed_macro_name, tissue_outline_save_end, tissue_outline);
    else
        [macro_processed_image, tissue_outline_im] = Step3_macrophage_post_processing_function(czi_image_folder, ...
            single_image_name, pre_processed_folder_name, post_processed_folder_name, ...
            vessel_seg_im, macro_seg_im, processed_macro_name, tissue_outline_save_end);
    end
    fprintf('step 3 complete\n')
%     toc
    %Step4
    
    tic
    [cell_array, labeled_macro_im] = Step4_label_macrophages_binary_func(czi_image_folder, ...
        single_image_name, post_processed_folder_name, ...
        macro_processed_image, '_labeled_macrophages_2023.tiff');
    
    fprintf('step 4 complete\n')
%     toc
    
%     tic
    [vessel_processed_im] =  Step5_vessel_post_processing_2022_func(czi_image_folder, ...
        single_image_name, post_processed_folder_name,vessel_seg_im, processed_ves_str);
    
    fprintf('step 5 complete\n')
%     toc
%     
%     tic
    Step6A_loop_ves_ntec_without_meanint(czi_image_folder, single_image_name,...
        pre_processed_folder_name, results_file_name_px, results_file_name_um, vessel_processed_im,...
        org_part_image, neighbourhood_dir_name, vessel_analysis_dir_name, ...
        metadata_file_name);
    fprintf('step 6 complete\n')
%     toc
%     
%     tic
    [thresholded_np_image, mask_np_image, inverted_tissue_image]= Step7_nanopaticle_post_processing_function(czi_image_folder, ...
        single_image_name, post_processed_folder_name, tissue_outline_im, ...
        org_part_image, nanoparticle_image_save_str, inverted_tissue_outline_name, thresholded_np_name );
    
    fprintf( 'step 7 complete\n')
%     toc
%   
    for rad_int = 1:size(radius_list,2)
        radius_um = radius_list(rad_int);
        STEP8_Identify_np_regions(czi_image_folder, single_image_name, ...
         post_processed_folder_name, np_region_dir_name, inverted_tissue_image, ...
        mask_np_image, radius_um, px_per_um )
    end


%     toc
    fprintf('completed_single_image')
end

toc