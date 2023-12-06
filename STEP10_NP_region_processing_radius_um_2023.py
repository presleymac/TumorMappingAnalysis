import re
import os
import time
import timeit

import numpy as np
import pandas as pd 
import skimage
from skimage import io
import tkinter as TK
from tkinter import filedialog #allows us to select the image of interest. 


def find_top_folder():
    '''find_top_folder
    Uses tkinter user interface to allow user to select the folder they want to
    process
    
    Returns:
    file_path (str): Path to the selected folder as a string.
    
    '''
    root = TK.Tk()
    root.fileName = filedialog.askdirectory()
    file_path = str(root.fileName)
    root.destroy()
    return(file_path)

def get_image_folder_list(folder):
    '''get_image_folder_list
    
    Arguments
    folder(str): Path to the folder in which to get a list of all subfolders
    
    Returns:
    image_folder_list(list): List of all subfolder in folder
    '''
    
    
    image_folder_list = next(os.walk(folder))[1]
    return(image_folder_list)

def read_image(image_path):
    ''' read_image
    
    Reads an image path into an np.ndarray.
    
    Arguments
    image_path(str): full path to image to be loaded
    
    Returns:
    image(np.ndarray): image as a numpy nd.array.
    
    '''
    
    image = io.imread(image_path)
    return(image)

def get_image_and_data_paths(folder, image_name, radius_um):    
    '''get_image_and_data_paths:
    
    Generates list of paths to the various images and files needed to run 
    analysis.
    
    Arguments:
    folder(str): Path to the top folder
    
    image_name(str):  String indicating the name of the image to generate sub images for
    
    radius_um:  Size of nanoparticle region radius to load paths for.
    
    Returns:
    np_image_path(str): path to processed nanoparticle image
    np_area_mask_path(str): path to nanoparticle region image
    macrophage_image_path(str) path to processed macrophage image
    tissue_image_path(str):Path to tissue image where tissue is given a value 1
    vessel_image_path(str):  path to processed vessel image
    vessel_data(str): path to analyzed vesel csv file generated in MATLAB
    
    '''

    np_image_path = (f'{folder}/{image_name}/post_processed_2023/'
                     f'{image_name}_post_processed_nanoparticles_2023.tiff')
    
    np_area_mask_path = (f'{folder}/{image_name}/Neighbourhood_Analysis_2023/NP_Region_Images/'
                     f'{image_name}_{radius_um}um_radius_np_regions_tissue_cropped_2023.tiff')
    
    macrophage_image_path = (f'{folder}/{image_name}/post_processed_2023/'
                     f'{image_name}_labeled_macrophages_2023.tiff')   
    
    tissue_image_path = (f'{folder}/{image_name}/post_processed_2023/'
                     f'{image_name}_tissue_outline_inverted_2023.tiff')     
    
    vessel_image_path =(f'{folder}/{image_name}/Neighbourhood_Analysis_2023/'
                        f'Vessel_Analysis_2023/{image_name}_labeled_vess_segs.tif')
    
    vessel_data = (f'{folder}/{image_name}/Neighbourhood_Analysis_2023/'
                    f'Vessel_Analysis_2023/{image_name}_complete_vessel_data_um_2023.csv')
  
    return(np_image_path, np_area_mask_path, macrophage_image_path, 
           tissue_image_path, vessel_image_path, vessel_data )
    
def load_data(folder, image_name, radius_um):
    '''load_data
    
    loads all the data necessary for analysis.
    
    Arguments:
    folder(str): Path to the top folder
    
    image_name(str):  String indicating the name of the image to generate sub images for
    
    radius_um:  Size of nanoparticle region radius to load paths for.
    
    Returns: 
    vessel_df(pd.DataFrame):  vessel data loaded as a pandas DataFrame
    np_im(np.ndarray): Nanoparticle image loaded into numpy nd.arrray
    np_mask_im(np.ndarray): Nanoparticle region postive mask image loaded into numpy nd.arrray 
    vessel_im(np.ndarray) Vessel image loaded into numpy nd.arrray
    macro_im (np.ndarray) Macrophage image loaded into numpy nd.arrray
    tissue_im(np.ndarray). Tissue image loaded into numpy nd.arrray where tissue = 1
    outside_np_mask(np.ndarray) Nanoparticle negative regin mask image loaded into numpy nd.arrray
    
    '''
    #Get image/data paths
    (np_image_path, 
     np_area_mask_path, 
     macro_image_path,
     tissue_image_path, 
     vessel_image_path, 
     vessel_data) = get_image_and_data_paths(folder, image_name, radius_um)
    
    #load vessel_feature data
    vessel_df=pd.read_csv(vessel_data, index_col=0)
    vessel_df['Image_name'] = image_name
    
    #Read in all the images we will be working with 
    
    np_im = read_image(np_image_path)
    np_mask_im = read_image(np_area_mask_path)
    vessel_im =  read_image(vessel_image_path)
    macro_im = read_image(macro_image_path)
    tissue_im = read_image(tissue_image_path)
    
    not_np_space = np_mask_im == 0 
    outside_np_mask = tissue_im*not_np_space

    
    return(vessel_df, np_im, np_mask_im, vessel_im, macro_im, tissue_im, 
           outside_np_mask)

def mask_data(feature_im, positive_mask, negative_mask=''):
    '''mask data
    Masks the feature image with the masks using element wise multiplication
    
    Arguments: 
    feature_image (np.ndarray): Image containing the features to be masked
    positive_mask(np.ndarray): Image containing the mask
    negative_mask(np.ndarray): Image containing the inverse mask
    
    Returns: 
    positve_features(np.ndarray): Masked image
    '''
    
    positive_feature = feature_im*positive_mask
    
    if not negative_mask=='':
        negative_feature = feature_im*negative_mask
        return(postive_feature), negative_feature
    else:
        return(positive_feature)

def quantify_features(folder, image_name, vessel_df, np_im, np_mask_im, vessel_im, macro_im, tissue_im, 
           outside_np_mask):
    
    ''' quantify_features
    
    Quantifies features inside nanoparticle postive and negative regions
    
    Arguments:
    folder (str): Path to image folder
    image_name (str): Name of image folder
    vessel_df(pd.DataFrame): Dataframe containing vessel analysis data
    np_im: (np.ndarray): Original post processed nanoparticle image
    np_mask_im (np.ndarray): Nanoparticle region image
    vessel_im(np.ndarray): Post processed vessel image
    macro_im(np.ndarray): Post processed macrophage image
    tissue_im(np.ndarray): Post processed tissue image
    outside_np_mask(np.ndarray): inverse of np_mask_image (Nanoparticle negative 
                                region mask)
    
    Returns:
    summary_df(pd.DataFrame): DataFrame containing a summary of all the features 
                                inside and outside of nanoparticle regions (full region 
                                sumary)
    vessel_outside_feature_df(pd.DataFrame): DataFrame containing all the vessel 
                                            segments classified as outside of nanoparticle 
                                            regions. (per vessel seg)
    
    
    vessel_inside_feature_df(pd.DataFrame): DataFrame containing all the vessel 
                                            segments classified as inside of nanoparticle 
                                            regions. (per vessel seg)
    
    
    
    '''

    
    #Use masks to isloate features in subset of tissue space
    vessel_inside_np = np_mask_im * vessel_im
    vessel_outside_np = outside_np_mask * vessel_im
    
    #Determine which vessels spread across np space boundaries
    shared_ves =np.intersect1d(np.unique(vessel_inside_np), np.unique(vessel_outside_np))
    
    #Find vessels inside and outside NP space and get a count of how many 
    # pixels are in each space
    vessel_inside_labels, vessel_inside_counts = np.unique(vessel_inside_np, return_counts=True)
    vessel_outside_labels, vessel_outside_counts = np.unique(vessel_outside_np, return_counts=True)

    #Generate a dataframe with the vessel label as index and the corresponding counts
    vessel_any_inside_df = pd.DataFrame({'label': vessel_inside_labels, 'inside_counts': vessel_inside_counts})
    vessel_any_outside_df = pd.DataFrame({'label': vessel_outside_labels, 'outside_counts': vessel_outside_counts})
    vessel_any_inside_df =  vessel_any_inside_df.set_index('label')
    vessel_any_outside_df =  vessel_any_outside_df.set_index('label')
    
    #Make a df containing all the vessels and the corresponding count inside and 
    # outside NP space. This allows for an easy comparison of label counts
    all_vessel = vessel_any_inside_df.merge(vessel_any_outside_df, how='outer', left_index=True, right_index=True)
    all_vessel = all_vessel.fillna(0)
    
    # Assign vessels to inside or outside np region. If 25% or more pixes are 
    # inside the np space we classify this as an "inside" feature. 
    # Otherwise this is an outside feature.
    inside_vessel_locations = (np.where(all_vessel['inside_counts']>= all_vessel['outside_counts']))
    outside_vessel_locations = np.where(all_vessel['inside_counts']< all_vessel['outside_counts'])

    # Find unique labels in the original image and corresponding pixel counts
    total_vessel_label, total_vessel_count = np.unique(vessel_im, return_counts=True)
    total_vessel_df =   pd.DataFrame({'label': total_vessel_label, 'counts': total_vessel_count})
    total_vessel_df = total_vessel_df.set_index('label')
    
    # Generate df of just vessels classified as inside or outside based in the 
    # 25% condition using the list generated above.

    vessel_inside_df = total_vessel_df.loc[inside_vessel_locations] 
    vessel_outside_df = total_vessel_df.loc[outside_vessel_locations]
    
    #try except to remove vessels labeled 0 as this is negative space 
    try:
        vessel_inside_df = vessel_inside_df.drop(index=0)
    except:
        pass
    
    try:
        vessel_outside_df = vessel_outside_df.drop(index=0)
    except:
        pass
    #Determine number of vessels in each category 
    total_inside_vessel_count =  len(vessel_inside_df)
    total_outside_vessel_count =  len(vessel_outside_df)

    #compare total pixel counts of vessels as classified by the 25% criteria to 
    #the number of pixels originally present inside the np space
    total_vessel_pixels_classified_inside = vessel_inside_df['counts'].sum()
    vessel_any_inside_df = vessel_any_inside_df.drop([0])
    total_vessel_pixels_inside = vessel_any_inside_df['inside_counts'].sum()
    
    total_vessel_pixels_classified_outside = vessel_outside_df['counts'].sum()
    vessel_any_outside_df = vessel_any_outside_df.drop([0])
    total_vessel_pixels_outside = vessel_any_outside_df['outside_counts'].sum()
    
    #Count the number of vessels in each category
    total_inside_vessel_count =  len(vessel_inside_df)
    total_outside_vessel_count =  len(vessel_outside_df)

    ### Repeat Steps for macrophages ###
    
    
    #Use masks to isloate features in subset of tissue space
    macro_inside_np = np_mask_im * macro_im
    macro_outside_np = outside_np_mask * macro_im
    
    #Determine which macrophages spread across np space boundaries
    shared_macro =np.intersect1d(np.unique(macro_inside_np), np.unique(macro_outside_np))
    #Find macro inside and outside NP space and get a count of how many 
    # pixels are in each space
    macro_inside_labels, macro_inside_counts = np.unique(macro_inside_np, return_counts=True)
    macro_outside_labels, macro_outside_counts = np.unique(macro_outside_np, return_counts=True)
    #Generate a dataframe with the macro label as index and the corresponding counts
    macro_any_inside_df = pd.DataFrame({'label': macro_inside_labels, 'inside_counts': macro_inside_counts})
    macro_any_outside_df = pd.DataFrame({'label': macro_outside_labels, 'outside_counts': macro_outside_counts})


    macro_any_inside_df =  macro_any_inside_df.set_index('label')
    macro_any_outside_df =  macro_any_outside_df.set_index('label')
    
    #Make a df containing all the pacrophagess and the corresponding count inside and 
    # outside NP space. This allows for an easy comparison of label counts    
    all_macro = macro_any_inside_df.merge(macro_any_outside_df, how='outer', left_index=True, right_index=True)
    all_macro = all_macro.fillna(0)
    
    # Assign vessels to inside or outside np region. If 25% or more pixes are 
    # inside the np space we classify this as an "inside" feature. 
    # Otherwise this is an outside feature.   
    inside_macro_locations = (np.where(all_macro['inside_counts']>= all_macro['outside_counts']))
    outside_macro_locations = np.where(all_macro['inside_counts']< all_macro['outside_counts'])

    # Find unique labels in the original image and corresponding pixel counts
    total_macro_label, total_macro_count = np.unique(macro_im, return_counts=True)
    total_macro_df =   pd.DataFrame({'label': total_macro_label, 'counts': total_macro_count})
    total_macro_df = total_macro_df.set_index('label')
    
    # Generate df of just vessels classified as inside or outside based in the 
    # 25% condition using the list generated above.
    macro_inside_df = total_macro_df.loc[inside_macro_locations] 
    macro_outside_df=total_macro_df.loc[outside_macro_locations]
    

    #try except to remove vessels labeled 0 as this is negative space     
    
    try:
        macro_inside_df = macro_inside_df.drop(index=0)
    except:
        pass
    
    try:
        macro_outside_df = macro_outside_df.drop(index=0)
    except:
        pass
    
    #Determine number of vessels in each category  
    total_inside_macro_count =  len(macro_inside_df)
    total_outside_macro_count =  len(macro_outside_df)
    #compare total pixel counts of vessels as classified by the 25% criteria to 
    #the number of pixels originally present inside the np space
    total_macro_pixels_classified_inside = macro_inside_df['counts'].sum()
    macro_any_inside_df = macro_any_inside_df.drop([0])
    total_macro_pixels_inside = macro_any_inside_df['inside_counts'].sum()
    #Count the number of vessels in each category
    total_macro_pixels_classified_outside = macro_outside_df['counts'].sum()
    macro_any_outside_df= macro_any_outside_df.drop([0])
    total_macro_pixels_outside = macro_any_outside_df['outside_counts'].sum()


    no_np_tissue = outside_np_mask
    
    total_np_empty_tissue_pix = no_np_tissue.sum()
    total_tissue_pix = tissue_im.sum()
    np_tissue_pix = np_mask_im.sum()
    
    total_inside_macro_num = len(macro_inside_df)
    total_outside_macro_num = len(macro_outside_df)
    
    ves_outside_list = np.unique(vessel_outside_df.index)
    ves_inside_list = np.unique(vessel_inside_df.index)
    
    vessel_outside_feature_df = vessel_df.loc[vessel_df['ves_label'].isin(ves_outside_list)].reset_index(drop=True)
    vessel_inside_feature_df =  vessel_df.loc[vessel_df['ves_label'].isin(ves_inside_list)].reset_index(drop=True)
    
    pct_inside_tissue_macro = total_macro_pixels_classified_inside/np_tissue_pix*100
    pct_outside_tissue_macro = total_macro_pixels_classified_outside/total_np_empty_tissue_pix*100
    
    pct_inside_tissue_vessel = total_vessel_pixels_classified_inside/np_tissue_pix*100
    pct_outside_tissue_vessel =  total_vessel_pixels_classified_outside/total_np_empty_tissue_pix*100
    
    total_np_intensity_inside = (np_im*np_mask_im).sum()
    total_np_intensity_outside = (np_im * outside_np_mask).sum()
    
    summary_df= pd.DataFrame({'np_positive_region': [0,1],
                              'image_name':[image_name,image_name],
                               'total_tissue_pixels':[total_np_empty_tissue_pix,
                                                      np_tissue_pix], 
                               'total_macro_classified_pix':[total_macro_pixels_classified_outside,
                                                             total_macro_pixels_classified_inside], 
                               'pre_classification_macro_pix': [total_macro_pixels_outside,total_macro_pixels_inside], 
                               'macro_count':[total_outside_macro_count,total_inside_macro_count], 
                               'total_ves_classified_pix': [total_vessel_pixels_classified_outside,total_vessel_pixels_classified_inside], 
                               'pre_classification_ves_pix': [total_vessel_pixels_outside,total_vessel_pixels_inside],  
                               'vessel_count': [total_outside_vessel_count,total_inside_vessel_count],
                               'pct_macro_in_tissue':[pct_outside_tissue_macro, pct_inside_tissue_macro],
                               'pct_ves_in_tissue':[pct_outside_tissue_vessel, pct_inside_tissue_vessel],#})
                               'np_intensity':[total_np_intensity_outside, total_np_intensity_inside] })


    return(summary_df, vessel_outside_feature_df, vessel_inside_feature_df)

def save_single_image_dfs(folder, image_name, summary_df, vessel_outside_feature_df, vessel_inside_feature_df, radius_um, units='um'):
    
    '''save_single_image_dfs
    
    saves dataframes for an individual image
    
    Arguments: 
    folder(str): Path to folder
    image_name(str): Image folder name
    summary_df(pd.DataFrame):  Full nanoparticle region summary values for image 
    vessel_outside_feature_df(pd.DataFrame): DataFrame containing all the vessel 
                                            segments classified as outside of nanoparticle 
                                            regions. (per vessel seg)
    vessel_inside_feature_df(pd.DataFrame): DataFrame containing all the vessel 
                                            segments classified as inside of nanoparticle 
                                            regions. (per vessel seg)
    radius_um(int): Radius size
    units='um'(str):  Optional argument which specifies the units of the radius 
                    input
    
    '''
    
    save_folder_path = f'{folder}/{image_name}/Neighbourhood_Analysis_2023/NP_Regions_2023'
    
    if not os.path.exists(save_folder_path):
        os.mkdir(save_folder_path)
    
    save_summary_path = f'{save_folder_path}/{image_name}_{radius_um}um_radius_np_region_summary_{units}.csv'
    save_outside_ves =  f'{save_folder_path}/{image_name}_{radius_um}um_radius_vessels_outside_np_regions_{units}.csv'
    save_inside_ves =  f'{save_folder_path}/{image_name}_{radius_um}um_radius_vessels_inside_np_regions_{units}.csv'
    
    summary_df.to_csv(save_summary_path)
    vessel_outside_feature_df.to_csv(save_outside_ves)
    vessel_inside_feature_df.to_csv(save_inside_ves)
    
def save_ves_summary_combo_df(folder, image_name, vessel_outside_feature_df, vessel_inside_feature_df,radius_um, units='um'):
    '''save_ves_sumary_combo_df
    Saves the combination DataFrames
    
    Arguments:
    folder(str): Path to folder
    image_name(str) Image folder name
    vessel_outside_feature_df(pd.DataFrame): DataFrame containing all the vessel 
                                            segments classified as outside of nanoparticle 
                                            regions with summary values. (per vessel seg)
    vessel_inside_feature_df(pd.DataFrame): DataFrame containing all the vessel 
                                            segments classified as inside of nanoparticle 
                                            regions with summary values. (per vessel seg)
    radius_um (int):  Radius size specified as an integer
    units='um'(str):  Optional argument which specifies the units of the radius 
                      input
    
    '''
    
    
    save_folder_path = f'{folder}/{image_name}/Neighbourhood_Analysis_2023/NP_Regions_2023'
    save_outside_ves =  f'{save_folder_path}/{image_name}_{radius_um}um_radius_vessels_outside_np_regions_with_summary_values_{units}.csv'
    save_inside_ves =  f'{save_folder_path}/{image_name}_{radius_um}um_radius_vessels_inside_np_regions_with_summary_values_{units}.csv'
    save_complete_ves = f'{save_folder_path}/{image_name}_{radius_um}um_radius_complete_vessels_with_summary_values_{units}.csv'
    
    full_image_summary_df =  pd.concat([vessel_outside_feature_df, vessel_inside_feature_df])
    full_image_summary_df.to_csv(save_complete_ves)
    vessel_outside_feature_df.to_csv(save_outside_ves)
    vessel_inside_feature_df.to_csv(save_inside_ves)
           
def prep_for_combo_df(summary_df, outside_ves_df, inside_ves_df):
    '''prep_for_comb_df
    
    Prepares data from quantify features function to be combined into one DataFrame
    
    Arguments:
    summary_df(pd.DataFrame): DataFrame containing a summary of all the features 
                                inside and outside of nanoparticle regions (full region 
                                sumary)
    vessel_outside_feature_df(pd.DataFrame): DataFrame containing all the vessel 
                                            segments classified as outside of nanoparticle 
                                            regions. (per vessel seg)
    
    
    vessel_inside_feature_df(pd.DataFrame): DataFrame containing all the vessel 
                                            segments classified as inside of nanoparticle 
                                            regions. (per vessel seg)
                                            
                                            
    Returns:
    outside_ves (pd.DataFrame): Vessels in Nanoparticle Negative regions with 
                                additional summary values added 
    inside_ves(pd.DataFrame): Vessels in Nanoparticle positive regions with 
                            additional summary values added 
    
    '''
    
    summary_to_add = summary_df.drop(columns = ['image_name'])
    column_names = summary_to_add.columns
    outside_summary = summary_to_add.iloc[0]
    inside_summary = summary_to_add.iloc[1]
    
    outside_ves= outside_ves_df.copy()
    inside_ves = inside_ves_df.copy()
    for column_int in range(len(column_names)):
        outside_ves[column_names[column_int]] = outside_summary[column_int] 
        inside_ves[column_names[column_int]] = inside_summary[column_int]
    
    return(outside_ves, inside_ves)

def save_full_folder_dataframe(full_folder_df,
                               top_folder_path,
                               radius_um, 
                               ending='um_radius_full_folder_vessel_data_um_2023.csv'):
    
    time_pt = folder_path.split('\\')[-1]
    save_folder =  f'{top_folder_path}/Summary_Data_2023'
    
    if not os.path.exists(save_folder):
        os.mkdir(save_folder)
    
    save_path =f'{save_folder}/{time_pt}_{radius}{ending}'
    full_folder_df.to_csv(save_path)

def subfolder_name_list(folder, included_list=['MSC', 'UT']):
    '''subfolder_name_list
    Generates a list of image subfolders in the top level folder that start 
    with the strings listed in the included list
    
    Arguments:
    folder (str): Path to the folder you want to generate a list of subfolders from.
    
    included_list(list): List of strings that indicate the start of names of the 
                        folders that should be included in the final subfolder list.
                        
    Returns:
    
    image_subfolder_list (list):  list of all the subfolders in a folder meeting 
                                    the naming condition.
    
    '''
    
    subfolder_name_list =  get_image_folder_list(folder)
    image_subfolder_list = []
    for image_string in included_list:
        for folder_name in subfolder_name_list:
            if folder_name.startswith(image_string):
                image_subfolder_list.append(folder_name)
    
    return(image_subfolder_list)
    

if __name__ == '__main__':

    #Specify the nanoparticle region radii to loop over,
    #Select the folder to analyze - this should be a time point folder (Read me folder level 3)
    
    radius_list = [25,50,100]
    folder =  find_top_folder()
    image_folder_list = subfolder_name_list(folder, included_list =['MSC','UT'])
    
    # Generate dataframe to store all the data in
            
    full_folder_df = pd.DataFrame()
    
    #Loop over image folders and print the name of the image being analyzed.
    for image_name in image_folder_list:
        print(image_name)

        #For each radius listed above load the images and ge

        for radius in radius_list:
            full_folder_df = pd.DataFrame()
            image_data_tuple = load_data(folder, image_name, radius)
            vessel_df = image_data_tuple[0] 
            np_im = image_data_tuple[1] 
            np_mask_im = image_data_tuple[2] 
            vessel_im = image_data_tuple[3] 
            macro_im = image_data_tuple[4] 
            tissue_im = image_data_tuple[5] 
            outside_np_mask= image_data_tuple[6] 
        
            single_image_summary, outside_ves, inside_ves = quantify_features(folder, 
                                                                            image_name, 
                                                                            vessel_df, 
                                                                            np_im, 
                                                                            np_mask_im, 
                                                                            vessel_im, 
                                                                            macro_im, 
                                                                            tissue_im, 
                                                                            outside_np_mask)
            
            outside_ves_sumamry_values, inside_ves_sumamry_values = prep_for_combo_df(single_image_summary, outside_ves, inside_ves)
        
            
            
            save_single_image_dfs(folder, image_name, single_image_summary, outside_ves, inside_ves,radius)
            
            save_ves_summary_combo_df(folder, image_name, outside_ves_sumamry_values, inside_ves_sumamry_values, radius)
