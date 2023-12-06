

mergestructs = @(x,y) cell2struct([struct2cell(x),struct2cell(y)],fieldnames(x),1);
matlab_folder = pwd;
image_seg_folder = uigetdir;
cd (image_seg_folder)

pre_processed_folder_name = 'pre_processing_2022';
post_processed_folder_name = 'post_processed_2023';


inverted_tissue_outline_name = '_tissue_outline_inverted_2023.tiff';
processed_macro_str = '_post_processed_macrophages_2023.tiff';
processed_ves_str = '_post_processed_vessels_2023.tiff';
nanoparticle_image_save_str = '_post_processed_nanoparticles_2023.tiff';


results_file_name_macro = '_macro_overlap_2023.xlsx';
results_file_name_vessel = '_vessel_overlap_2023.xlsx';

neighbourhood_dir_name = 'Neighbourhood_Analysis_2023';
overlap_analysis_dir_name = 'Overlap_Coefficients';



MSC = dir('MSC*');
UT = dir('UT*');
files = mergestructs(MSC,UT);

for image = 1:size(files,1)
    [~,single_image_name] = fileparts(files(image).name);
    
    image_folder_dir = strcat(image_seg_folder,'\', single_image_name);
    
    cd(image_folder_dir)
    tic
    %get paths and load images
    post_processed_dir =  strcat(image_folder_dir, '\', post_processed_folder_name);
    save_dir = strcat(image_folder_dir, '\', neighbourhood_dir_name, '\', overlap_analysis_dir_name);
    
    
    vessel_name = strcat(single_image_name, processed_ves_str);
    macro_name = strcat(single_image_name, processed_macro_str);
    tissue_outline_name = strcat(single_image_name, inverted_tissue_outline_name);
    np_image_name = strcat(single_image_name,nanoparticle_image_save_str); 
    
    cd(post_processed_dir)
    macro_im = imreadfast(macro_name);
    vessel_im = imreadfast(vessel_name);
    particle_im = imreadfast(np_image_name);
    
   if exist(save_dir, 'dir')~=7
        mkdir(save_dir)
   end 
   
   
    [pearson_coef_coloc, manders_overlap, k1, k2, M1, M2, mask_im, ...
    particle_threshold_image, bkg_corrected_particle_im] = ...
    tumour_mapping_correlation_coeff_fun(macro_im, particle_im);

    table_macro = table(pearson_coef_coloc, manders_overlap, k1, k2, M1, M2, mask_im, ...
    particle_threshold_image, bkg_corrected_particle_im)
    
    [pearson_coef_coloc, manders_overlap, k1, k2, M1, M2, mask_im, ...
    particle_threshold_image, bkg_corrected_particle_im] = ...
    tumour_mapping_correlation_coeff_fun(vessel_im, particle_im);

    table_vessel = table(pearson_coef_coloc, manders_overlap, k1, k2, M1, M2, mask_im, ...
    particle_threshold_image, bkg_corrected_particle_im)
    
    [pearson_coef_coloc, manders_overlap, k1, k2, M1, M2, mask_im, ...
    particle_threshold_image, bkg_corrected_particle_im] = ...
    tumour_mapping_correlation_coeff_fun(macro_im, vessel_im);
    
    table_macro_ves = table(pearson_coef_coloc, manders_overlap, k1, k2, M1, M2, mask_im, ...
    particle_threshold_image, bkg_corrected_particle_im)

end
    

