#process vessels only
import sys

import seaborn as sns
sns.set()
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_samples, silhouette_score, calinski_harabasz_score

from scipy.spatial import ConvexHull, Delaunay
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import kneed 
from kneed import KneeLocator
import tkinter as TK
from tkinter import filedialog #allows us to select the image of interest. 
import time
import timeit

import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np 
import re
import os
from skimage import io


def find_top_folder():
    """ find_top_folder
    
    This function can be used to open a file dialog box where the user can 
    select the folder the want the path for 
    
    Returns: 
    file_path (str): Returns the path to the folder as a string """
    
    root = TK.Tk()
    root.fileName = filedialog.askdirectory()
    file_path = str(root.fileName)
    root.destroy()
    return(file_path)
   

def vessel_data_pre_processing (ves_data_df,
                                ves_seg_analysis_df, 
                                image_name,
                                top_folder,
                                save_path=''):
    
   
   
   
    neigh_analysis = 'Neighbourhood_Analysis_2023'
    dataframe = 'Dataframes'
    ves_analysis = 'Vessel_Analysis_2023'
    
    
    
    if save_path == '':
        # print('Select the top folder')
        top_path = top_folder
        ves_save_path = f'{top_path}/{image_name}/{neigh_analysis}/{ves_analysis}'
        dataframe_save_path = f'{top_path}/{image_name}/{neigh_analysis}/{dataframe}'
    

    complete_vessel_df = pd.concat([ves_data_df, ves_seg_analysis_df], axis=1) 
    save_complete_vessel_data(complete_vessel_df, ves_save_path, image_name)

    return( complete_vessel_df)

def load_all_ves_data(top_folder,
                      single_image_name, 
                      results_ending='Vessel_analysis_results_um_2023.xlsx',
                      assignment_ending = '_correct_vessel_cluster_assignment.csv',
                      complete_ves_df = '_complete_vessel_data_2023.csv', 
                      vessel_summary = '_vessel_analysis_summary.csv',
                      ves_ending = '_vessel_percent_df.csv',
                      vessel_im_ending = '_labeled_vess_segs.tif'
                       ):
    
    '''load_all_ves_data
    Loads data from MATLAB output and adds a vessel label to the DataFrame and saves the results
    
    Arguments:
    top_folder(str): Path to time point folder
    single_image_name(str): Image folder name
    
    results_ending(str): Default value:'Vessel_analysis_results_um_2023.xlsx' specifies ending of the 
                    file name as output by MATLAB scripts
                    
    assignment_ending(str): Default value: '_correct_vessel_cluster_assignment.csv' ending of the path for the 
                            Dataframe being saved
    complete_ves_df(str): Default '_complete_vessel_data_2023.csv' ending of the path for the 
                            Dataframe being saved
    vessel_summary(str): Default '_vessel_analysis_summary.csv' ending of the path for the 
                            Dataframe being saved
    vessel_im_ending(str) Default '_labeled_vess_segs.tif' ending of the labelled vessel image being
                            loaded as an np.ndarray
    

    Returns: 
    
    ves_analysis_df (pd.DataFrame): Vessel data directly from MATLAB
    
    complte_ves_df(pd.DataFrame): Vessel data with labels    
    '''
    
    
    
    vessel_data_folder = 'Neighbourhood_Analysis_2023/Vessel_Analysis_2023'
    vessel_data_path = f'{top_folder}/{single_image_name}/{vessel_data_folder}'
        
    vessel_analysis_path = f'{vessel_data_path}/{single_image_name}{results_ending}'
    
    complete_ves_path = f'{vessel_data_path}/{single_image_name}{complete_ves_df}'
    vessel_summary_path = f'{vessel_data_path}/{single_image_name}{vessel_summary}'
    
    
    ves_analysis_df = pd.read_excel(vessel_analysis_path, engine='openpyxl')
    try:
        vessel_cluster_assignment = pd.read_csv(vessel_cluster_assin_path)
        complete_ves_df =  pd.read_csv(complete_ves_path)
        vessel_summary_df = pd.read_csv(vessel_summary_path)
    except:

        save_path = f'{vessel_data_path}/{single_image_name}'
        ves_im_path = f'{vessel_data_path}/{single_image_name}{vessel_im_ending}'
        ves_im = io.imread(ves_im_path)
        ves_labels = range(len(ves_analysis_df))
        ves_df = pd.DataFrame()
        ves_df['ves_label'] = ves_labels
        complete_ves_df =  vessel_data_pre_processing(ves_df,ves_analysis_df, single_image_name, top_folder)

        save_complete_vessel_data(complete_ves_df, 
                                  vessel_data_path, 
                                  single_image_name, 
                                  ending='_complete_vessel_data_um_2023.csv')
    
    return(ves_analysis_df, complete_ves_df)

def save_complete_vessel_data(complete_ves_df, 
                              save_path, 
                              image_name, 
                              ending='_complete_vessel_data_um_2023.csv'):
    
    '''save_complete_vessel_data
    Saves complete_ves_df as pd.DataFrame
    
    '''
    complete_ves_df.to_csv(f'{save_path}/{image_name}{ending}')
    
def save_vessel_summary_data(vessel_summary_df, 
                             save_path, 
                             image_name,  
                             ending='_vessel_analysis_summary_um_2023.csv'):
    
    '''save_vessel_summary_data
    
    Saves vessel data as .csv
    
    '''
    
    vessel_summary_df.to_csv(f'{save_path}/{image_name}{ending}')

def get_image_folder_list(folder):
    '''get_image_folder_list
    Loops over time point folder to generate a list of image subfolders
    
    Arguments: 
    folder(str): Path to the folder being analyzed 

    Returns: 
    imagefolder_list(list): List of image sumfolder names
    '''
    
    image_folder_list = next(os.walk(folder))[1]
    return(image_folder_list)

def subfolder_name_list(folder, included_list=['MSC', 'UT']):
    '''subfolder_name_list
    Loops over time point folder to generate a list of image subfolders that start with
    the names in included list.
    
    Arguments: 
    folder(str): Path to the folder being analyzed 
    included_lsit(list):  List of strings indicating image folder names to include in subfolder list
    
    Returns: 
    image_subfolder_list(list): List of image sumfolder names
    
    '''
    
    subfolder_name_list =  get_image_folder_list(folder)
    image_subfolder_list = []
    for image_string in included_list:
        for folder_name in subfolder_name_list:
            if folder_name.startswith(image_string):
                image_subfolder_list.append(folder_name)
    
    return(image_subfolder_list)
  
if __name__ == '__main__':
    
    folder_path  = find_top_folder()
    image_name = subfolder_name_list(folder_path)
    for image in image_name:
        print(image)
        ves_analysis_df, complete_ves_df = load_all_ves_data(folder_path, image)
