import re
import os
import time
import timeit

import numpy as np
import pandas as pd 
import skimage
from skimage import io
import tkinter as TK
from tkinter import filedialog 
import neighbourhood_analysis_functions as nf
from NP_region_processing_radius_um_2023 import subfolder_name_list



def merge_and_image_time_folder_summary_values(time_point_folder_path, radius):
    
    image_subfolder_list = subfolder_name_list(time_point_folder_path)
    time_pt_str =  time_point_folder_path.split('/')[-1]
    time_pt_float = float(time_pt_str.split('h')[0])
    part_path =time_point_folder_path.split('/')[1]
    part_size_name = part_path.split('\\')[-1]
    part_size_str = part_path.split('GNP')[-1]
    part_size_float = float(part_size_str.split('nm')[0])
    

    save_folder = (f'{time_point_folder_path}/'
                 f'Neighbourhood_Analysis_2023/Summary_Dataframes/')
                 
    if not os.path.exists(save_folder):
        os.mkdir(save_folder)
        
    save_path = (f'{save_folder}/{part_size_name}_{time_pt_str}'
                f'_{radius}um_radius_full_folder_sumamry_values.csv')             
    
    time_pt_summary_region_df = pd.DataFrame()
    for image in image_subfolder_list:
                    
        summary_region_path = (f'{time_point_folder_path}/{image}/'
                                f'Neighbourhood_Analysis_2023/NP_Regions_2023/'
                                f'{image}_{radius}um_radius_np_region_summary_um.csv')
 
        single_img_summary_regions = pd.read_csv(summary_region_path, index_col=[0])
        time_pt_summary_region_df = pd.concat([time_pt_summary_region_df, 
                                               single_img_summary_regions])
        
    time_pt_summary_region_df['time_pt'] = time_pt_float
    time_pt_summary_region_df['np_size'] = part_size_float
    time_pt_summary_region_df['radius'] = radius
    
    time_pt_summary_region_df.to_csv(save_path)
    
    

def merge_and_save_np_size_region_summary(size_folder_path, radius):
    time_pt_list =  next(os.walk(size_folder_path))[1]
    time_pt_list = [x for x in time_pt_list if x.endswith('h')]
    size = size_folder_path.split('/')[-1]
    
    save_folder = (f'{size_folder_path}/Neighbourhood_Analysis_2023/'
                   'Summary_Dataframes')
    
    if not os.path.exists(save_folder):
        os.mkdir(save_folder)              
    
    save_path = f'{save_folder}/{size}_{radius}um_radius_all_time_pts_summary_values.csv'
            
    size_df = pd.DataFrame()
    for time_pt in time_pt_list:
        time_folder_path = f'{size_folder_path}/{time_pt}'
   
        folder_summary = (f'{time_folder_path}/Neighbourhood_Analysis_2023/'
                          f'Summary_Dataframes/{size}_{time_pt}_{radius}um'
                          '_radius_full_folder_sumamry_values.csv')

        time_pt_df = pd.read_csv(folder_summary, index_col=[0])
        size_df = pd.concat([size_df, time_pt_df])
    
    size_df.to_csv(save_path)    

def full_tumour_summary_dataframe(top_folder, radius):    
    
    save_folder = f'{top_folder}/Neighbourhood_Analysis_2023/Summary_Dataframe'
    if not os.path.exists(save_folder):
        os.mkdir(save_folder)
    
    save_name = f'{save_folder}/{radius}um_radius_full_tumour_summary_values.csv'    
    
    tumour_df = pd.DataFrame()
    size_folder_name_list = next(os.walk(top_folder_path))[1]
    size_folder_name_list = [x for x in size_folder_name_list if x.startswith('U87')]
    for size in size_folder_name_list:
        size_folder_path = f'{top_folder_path}/{size}'
        single_np_size = (f'{size_folder_path}/Neighbourhood_Analysis_2023/'
                          f'Summary_Dataframes/{size}_{radius}um_radius_all_time_pts_summary_values.csv')
        single_size_df = pd.read_csv(single_np_size, index_col=[0])
        
        tumour_df = pd.concat([tumour_df, single_size_df])
    
    tumour_df.to_csv(save_name)


if __name__ == '__main__':
    
    top_folder_path = nf.get_folder_path()
    size_folder_name_list = next(os.walk(top_folder_path))[1]
    size_folder_name_list = [x for x in size_folder_name_list if x.startswith('U87')]
    for size in size_folder_name_list:

        size_folder_path = f'{top_folder_path}/{size}'
        
        time_pt_list =  next(os.walk(size_folder_path))[1]
        time_pt_list = [x for x in time_pt_list if x.endswith('h')]
        radius_list = [25, 50, 100]
        time_point_path_list = []
        for time_pt in time_pt_list:
            time_folder_path = f'{size_folder_path}/{time_pt}'   
            time_point_path_list.append(time_folder_path)
            
    
        for time_pt_path in time_point_path_list:
            print(time_pt_path)
            for radius in radius_list:   
                print(radius) 
                merge_and_image_time_folder_summary_values(time_pt_path, radius) 
        
        
        for radius in radius_list:
            print(radius)
            merge_and_save_np_size_region_summary(size_folder_path, radius)
            full_tumour_summary_dataframe(top_folder_path, radius)


