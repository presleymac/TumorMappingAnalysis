function [vessel_processed_im] =  Step5_vessel_post_processing_2022_func(image_folder_path, ...
    image_name, post_processed_dir_name, ...
    vessel_seg_img, post_processed_vessel_name)
% Saves the path the to matlab folder (folder containing the script)
% User then selects the folder containing the image segmentations and
% the name of the images to be anaylyzed (UT or MSC)

image_seg_folder = image_folder_path;
cd (image_seg_folder)

%Section 1: Define the image folders and generate paths or folders as
%needed.

tic

% first image shortnames are identified as the name of each image
% subfolder in the original directory. A message is displayed to the
% user indicating which of the images is currently being processed we
% change directory to that image subfolder. Next subfolders where the
% simple segementation images are contained are identified.

display (['Processing Vessels ' image_name])
image_folder_dir = strcat(image_seg_folder,'\', image_name);
cd(image_folder_dir)

post_processing_dir=(strcat(image_folder_dir,'\', post_processed_dir_name));


%%
%Section 2: find pre_processed_images and begin loading the variables.
%Note tissue outlines are defined in the vessel images thus we load the
%macrophage and vessel images into the script.

%all vessels are stained with Cy3 = Channel 3
%all macrophages are stained with Cy5 = Channel 4
%Dapi was not used to label the tissue as not all tissues were DAPI
%stained.

cd(post_processing_dir)
vessel_processed_name = strcat(image_name, post_processed_vessel_name);
vessel_channel = vessel_seg_img;
ves_thresh_bin = (vessel_channel == 1);
% se2 = strel('sphere', 2);
% remove_sm_noise = bwareaopen(vessels, 26, 26);

% dil_noise_free = imdilate(remove_sm_noise, se2);
% er = erosion( dil_noise_free, 3, 'elliptic');
% ves_processed = erosion(er, 3, 'elliptic');

% ves_processed_uint16 = uint16(ves_processed);


se3 = strel('sphere',3);
ves_thresh_close = imclose(ves_thresh_bin,se3);
ves_processed_uint16 = uint16(ves_thresh_close);


vessel_file_name= vessel_processed_name;
clear options;
options.overwrite = true;
options.compress = 'lzw';
saveastiff(uint16(ves_processed_uint16), vessel_file_name, options);

vessel_processed_im = ves_processed_uint16;
toc
end

