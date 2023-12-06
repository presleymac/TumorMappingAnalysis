function [thresholded_np_image, mask_np_image, inverted_tissue_image]= Step7_nanopaticle_post_processing_function(image_folder_path, ...
    image_name, post_processed_dir_name, tissue_outline_im, ...
    original_nanoparticle_image, nanoparticle_image_save_str, inverted_tissue_outline_name, thresholded_np_name )

% Saves the path the to matlab folder (folder containing the script)
% User then selects the folder containing the image segmentations and
% the name of the images to be anaylyzed (UT or MSC)


image_seg_folder = image_folder_path;
cd (image_seg_folder)

% first image shortnames are identified as the name of each image
% subfolder in the original directory. A message is displayed to the
% user indicating which of the images is currently being processed we
% change directory to that image subfolder. Next subfolders where the
% simple segementation images are contained are identified.

display (['Processing Nanoparticles ' image_name])
image_folder_dir = strcat(image_seg_folder,'\', image_name);
cd(image_folder_dir)

post_processing_dir=(strcat(image_folder_dir,'\', post_processed_dir_name));


np_processed_name = strcat(image_name, nanoparticle_image_save_str);

tissue_inverted =  strcat(image_name, inverted_tissue_outline_name);

thresholded_np_save_name =  strcat(image_name, thresholded_np_name);

cd(image_folder_dir);
NP_im = (original_nanoparticle_image);
cd(post_processing_dir)
thresh_val = graythresh(NP_im);
thresholded_np = imbinarize(NP_im, thresh_val*2);
thresh_np = uint16(thresholded_np);
tissue_im = tissue_outline_im;
invert_tissue = tissue_im == 0;
invert_tissue = uint16(invert_tissue);
mask_np = thresh_np.*NP_im;




tissue_inverted_name = tissue_inverted;
clear options;
options.overwrite = true;
options.compress = 'lzw';
saveastiff(uint16(invert_tissue), tissue_inverted_name, options);



mask_name = np_processed_name;
clear options;
options.overwrite = true;
options.compress = 'lzw';
saveastiff(uint16(mask_np), mask_name, options);

threshold_name = thresholded_np_save_name;
clear options;
options.overwrite = true;
options.compress = 'lzw';
saveastiff(uint16(thresholded_np), threshold_name, options);

thresholded_np_image = thresholded_np;
mask_np_image =  mask_np;
inverted_tissue_image = invert_tissue;


end

