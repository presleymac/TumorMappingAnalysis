import os
import numpy as np
import pandas as pd
import raster_geometry as rg
import time
from skimage import io
import neighbourhood_analysis_functions as nf

def sphere_sampling(path_to_image_folder,
                    image_name,
                    tissue_im,
                    vessel_im,
                    macro_im,
                    np_im,
                    radius_um,
                    sample_number,
                    metadata_pix_size):


    '''sphere_sampling
    Randomly sample image to genereate list of center points at which spheres of 
    desired radius ize are generated. 
    
    
    Arguments: 
    path_to_image_folder (str): Path to the image folder
    image_name (str): Iamge folder name
    tissue_im (np.ndarray):  tissue image to sample
    vessel_im (np.ndarray): Vessel image to sample
    macro_im (np.ndarray): Macrophage image to sample
    np_im (np.ndarray): Nanoparticle image to sample
    radius_um(int): desired sphere radius in um
    sample_number(int):  number of randomly generated spheres to make
    metadata_pix_size(int): Relationship between image pixel size and um
    
    
    Returns:
    iod_df(pd.DataFrame):  DataFrame containing quantificatin of features in each
                            spherical sample. DataFrame contains the data used to 
                            calculate IOD in a later script.

    center_list(pd.DataFrame): DataFrame of sphere center points
    
    
    dataframe_save_path(str): path to the save folder where the dataframes are 
                                being saved
                    
    '''
    
    #dataframe to save results
    iod_df = pd.DataFrame()



    save_folder_path = f'{path_to_image_folder}/Neighbourhood_Analysis_2023/Summary_Dataframes'

    if not os.path.exists(save_folder_path):
        os.mkdir(save_folder_path)

    dataframe_save_path = f'{save_folder_path}/{image_name}_{sample_number}x_{radius_um}um_sphere_sampling_index_of_dispersion.csv'
    center_pt_path = f'{save_folder_path}/{image_name}_{sample_number}x_{radius_um}um_sphere_sampling_center_pts.csv'


    radius = radius_um/metadata_pix_size
    center_list = []
    z_size = tissue_im.shape[0]
    y_size = tissue_im.shape[1]
    x_size = tissue_im.shape[2]

    image_shape = (z_size, y_size, x_size)
    mask_size_list = []
    tissue_pix_list = []
    macro_pix_list = []
    vessel_pix_list = []
    np_pix_list = []
    pct_tissue_list = []
    pct_macro_in_tis_list = []
    pct_ves_in_tis_list = []
    pct_np_in_tis_list = []
    z_list = []
    y_list = []
    x_list = []

    count = 0
    while count < sample_number:
        center = (
    np.random.randint(high=z_size, low=0),
    np.random.randint(high=y_size, low=0),
    np.random.randint(high=x_size, low=0)
    )

        center_list.append(center)
        count = count + 1
    center_df = pd.DataFrame(center_list, columns=['z', 'y', 'x'])
    center_df.to_csv(center_pt_path)
    for center_pt in center_list:

        # start_time = time.time()
        z = center_pt[0]
        y = center_pt[1]
        x = center_pt[2]
        sphere_mask = (rg.nd_superellipsoid(
                image_shape, radius, 2.0, center_pt, 3,
                rel_position=False, rel_sizes=False, smoothing=False)).astype(int)


        mask_pix = np.count_nonzero(sphere_mask)
        tissue_pix =np.count_nonzero(sphere_mask *tissue_im)
        if tissue_pix==0:
            continue
        else:
            macro_pix = np.count_nonzero(sphere_mask * macro_im)
            vessel_pix = np.count_nonzero(sphere_mask * vessel_im)
            np_pix = np.count_nonzero(sphere_mask * np_im)

        # print("--- %s seconds ---" % (time.time() - start_time))


        pct_tissue = tissue_pix/mask_pix*100


        pct_macro_in_tissue = macro_pix/tissue_pix*100
        pct_ves_in_tissue = vessel_pix/tissue_pix*100
        pct_np_in_tissue = np_pix/ tissue_pix*100

        z_list.append(z)
        y_list.append(y)
        x_list.append(x)
        mask_size_list.append(mask_pix)
        tissue_pix_list.append(tissue_pix)
        macro_pix_list.append(macro_pix)
        vessel_pix_list.append(vessel_pix)
        np_pix_list.append(np_pix)
        pct_tissue_list.append(pct_tissue)
        pct_macro_in_tis_list.append(pct_macro_in_tissue)
        pct_ves_in_tis_list.append(pct_ves_in_tissue)
        pct_np_in_tis_list.append(pct_np_in_tissue)

        # print("--- %s seconds ---" % (time.time() - start_time))
    iod_df['z'] = z_list
    iod_df['y'] = y_list
    iod_df['x'] = x_list
    iod_df['mask_size'] = mask_size_list
    iod_df['tissue_pixels'] = tissue_pix_list
    iod_df['macro_pixels'] = macro_pix_list
    iod_df['vessel_pixels'] = vessel_pix_list
    iod_df['np_pixels'] = np_pix_list
    iod_df['pct_tissue'] = pct_tissue_list
    iod_df['pct_macro_per_tissue'] = pct_macro_in_tis_list
    iod_df['pct_ves_per_tissue'] = pct_ves_in_tis_list
    iod_df['pct_np_per_tissue'] = pct_np_in_tis_list



    iod_df.to_csv(dataframe_save_path)

    return(iod_df, center_list, dataframe_save_path)


if __name__ == '__main__':

    #Select time point folder and generates a list of image folder names
    top_folder_path = nf.get_folder_path()
    subfolders = nf.generate_list_subfolders(top_folder_path)

    image_list = []

    for im in subfolders:
        if 'MSC' in im:
            image_list.append(im)
        # elif 'MSC162' in im:
        #     image_list.append(im)
        elif 'UT' in im:
            image_list.append(im)

    #specify sample number an the desired radius
    sample_num = 1000
    radius_um = int(input("enter radius:"))
    for image_folder_path in image_list:

        start_time = time.time()
        image_name = image_folder_path.split('\\')[-1]

        print('Sampling ',image_name, 'and saving sample df ')
        #find paths and load data
        processed_folder = f'{image_folder_path}/post_processed_2023'

        tissue_path = f'{processed_folder}/{image_name}_tissue_outline_inverted_2023.tiff'
        vessel_path = f'{processed_folder}/{image_name}_post_processed_vessels_2023.tiff'
        macro_path = f'{processed_folder}/{image_name}_post_processed_macrophages_2023.tiff'
        np_path = f'{processed_folder}/{image_name}_post_processed_nanoparticles_2023.tiff'
        metadata_path = f'{image_folder_path}/{image_name}_iso_info.csv'
        meta_data = pd.read_csv(metadata_path)
        pixel_size = round(float(meta_data['newphys'][0])*(10**6))

        tissue_im = io.imread(tissue_path)
        vessel_im = io.imread(vessel_path)
        macro_im = io.imread(macro_path)
        np_im = io.imread(np_path)


        #Generates IOD_df, center_points and path to data
        df, pts, pth =  sphere_sampling(image_folder_path,
                                        image_name,
                                        tissue_im,
                                        vessel_im,
                                        macro_im,
                                        np_im,
                                        radius_um,
                                        sample_num,
                                        pixel_size)
        print("--- %s seconds ---" % (time.time() - start_time))
