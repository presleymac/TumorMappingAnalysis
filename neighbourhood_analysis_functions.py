
#neighbourhood analysis functions



import os
import numpy as np
import pandas as pd
import itertools
import collections
import matplotlib.pyplot as plt
from skimage import io 
from statsmodels.stats import weightstats
import skgstat
from skimage.filters import threshold_otsu
from statsmodels.stats import weightstats
from scipy import stats
from scipy.stats import zscore
from scipy.ndimage import distance_transform_edt as transform_edt
import tkinter as TK

from tkinter import filedialog
import time
from copy import deepcopy


def get_folder_path(print_statement=False):
    """ get_folder_path
    
    This function can be used to open a file dialog box where the user can 
    select the folder the want the path for 
    
    Arguments:
    print_statement (logical): If true print statments are generated which 
                                print the folder path to the terminal
    
    Returns: 
    path (str): Returns the path to the folder as a string """
    
    root = TK.Tk()
    root.fileName = filedialog.askdirectory(title = ("""Select the folder to 
                                                     generate the path for"""))
    directory = str(root.fileName)
    top_level_folder = directory.replace('/', '\\')
    root.destroy()
    path = top_level_folder
    if print_statement==True:
        print(f'Image folder path is {path}\n')
    return(path)
        
        
def post_processed_folder_path(animal_folder_path="",
                                post_processed_folder_name="post_processed_2022",
                                print_statment=False):
    """ post_processed_folder_path
    
    This function generates the path to the post_processed_folder. This folder 
    should contain the Macrophage, Vessel and Nanoparticle images in their post
    processed (ready to analyze form). This function takes only optional 
    arguements. The animal_folder_path should be given or the user will be prompted 
    to select it. This function is designed to easily be incorportated into loops
    allowing the user to loop over all the animal folders in for a given time 
    point. 
    
    See also loop_all_folders()
    
    Arguments: 
    animal_folder_path (str): Default = "". This argument is the path to the 
                            image_folder being analyzed. If no path is specifed 
                            the user is prompted to select one. Note this should
                            be of the form "C:\\loction_of_data\\tumour_mapping
                            \\timepoint\\ANIMAL_NAME"
    
    post_processed_folder_name (str): Default value = "post_processed_2022".
                                    This should specify the name of the folder
                                    containing post processed Macrophage, 
                                    Vessel and Nanoparticle images. If another 
                                    folder name is not specified the default is 
                                    assumed.

    print_statement (logical): Default=False. If true print statments are
                            generated which print the folder path to the 
                            terminal.
                            
    Returns:
    post_processed_path (str): Returns the path to the post processed folder.

    """
    
    if animal_folder_path == "":
        animal_folder_path = get_folder_path()
    else:
        animal_folder_path = animal_folder_path
    
    post_processed_path = f'{animal_folder_path}\\{post_processed_folder_name}'
    
    if print_statment==True:
        print(f'Post processed image path {post_processed_path}\n')
        
    return(post_processed_path)

def tissue_outline_folder_path(animal_folder_path="", 
                          tissue_outline_folder_name='post_processed_2022',
                          print_statement=False):
    
    
    """ tissue_outline_folder_path:
    
    This function generates the path to the folder containing the tissueoutline.
    A default folder name is used if another is not specified by the user.
    
    Arguments:
    animal_folder_path (str): Default = "". This argument specifies the path to 
                            the image_folder being analyzed. If no path is specifed 
                            the user is prompted to select one. Note this should
                            be of the form "C:\\loction_of_data\\tumour_mapping
                            \\timepoint\\ANIMAL_NAME"
                            
    tissue_outline_folder_path (str): Default = 'tissueoutline'. This argument
                                    is used to specify the name of the 
                                    folder containing the tissue outline.
    
    print_statement(logical): Default=False. If true print statments are
                            generated which print the folder path to the 
                            terminal.
                            
    Returns:
    tissue_outline_folder_path (str): The path to the folder containing the 
                                    tissue outline.

                            
    """
    
    if animal_folder_path == "":
        animal_folder_path = get_folder_path()
    else:
        animal_folder_path = animal_folder_path
         
    tissue_outline_folder_path =f'{animal_folder_path}\\{tissue_outline_folder_name}'
    
    if print_statement==True:
        print(f'Tissue outline folder path {tissue_outline_folder_path}\n')
    return(tissue_outline_folder_path)
    
def generate_list_subfolders(top_level_folder_path):
    
    if top_level_folder_path == "":
        top_path = get_folder_path()  
    else:
        top_path = top_level_folder_path

    paths = []
    im_stack_name_subfolder_list = next(os.walk(top_path))[1] 
    
    for folder in im_stack_name_subfolder_list:
        if folder.startswith('_'):
            continue
        else:
            folder_name= folder

            
            path = f'{top_path}\\{folder_name}'
            paths.append(path)

    return(paths)

def return_stack_name(path="",
                      print_statement=False):
    
    '''This function identifies the name of the tumour image being used (i.e. 
        MSC160-T-stack3) and the path to the master folder containing all images 
        for that tumour.
        
        Note: To use this in a loop feed the individual path into this function.
        
        Returns: 
        image_shortname_name(str): The name of the image (i.e MSC160-T-stack3)
        image_folder_path(str): The path to the master image folder
        '''
        #locate the image folder 
        
    if path == "":
        path = get_folder_path()
    else:
        path = path     
    
    image_shortname = path.split('\\')[-1]
    
    if print_statement==True:
        print(f'image shortname is:  {image_shortname}')
    
    return (image_shortname)    
    
def simple_path_gen (top_level_folder,
                    stack_name="",
                    macro_ending='_labeled_macrophages.tiff',
                    vessel_ending='_vessels_post_processed_2023.tiff',
                    np_ending='_particles_post_processed_2023.tiff',
                    tissue_outline_ending='_tissue_outline_2023.tiff'): 
    
    if top_level_folder == "":
        top_path = get_folder_path()
    else:
        top_path = top_level_folder
    
    if stack_name == "":
        stack_name = return_stack_name()
        
    animal_folder_path =  f'{top_path}\\{stack_name}'    
        
    post_processed_folder = post_processed_folder_path(animal_folder_path)
    tissue_outline_folder = tissue_outline_folder_path(animal_folder_path)
    image_shortname = return_stack_name(animal_folder_path)

    
    macro_image_path =f'{post_processed_folder}\\{image_shortname}{macro_ending}'
  
 
    vessel_image_path = f'{post_processed_folder}\\{image_shortname}{vessel_ending}'


    np_image_path = f'{post_processed_folder}\\{image_shortname}{np_ending}'

        
   
    tissue_outline_path = f'{tissue_outline_folder}\\{image_shortname}{tissue_outline_ending}'

    return( macro_image_path, vessel_image_path, np_image_path, tissue_outline_path)  
      
def loop_all_folders(top_level_folder_path="",
                        animal_name=""):
    
    ''' loop_all_folders
    
    This folder will loop over all animal folders in the top_level_folder_path. 
    If an animal_name is specified only one folder will be processed. The result
    id a list of all image sub folders (I.e inside of a time folder). 
    
    Arguments
    
    top_level_folder_path (str): This specifies the top level folder. This should
                                be the path to the folder containing all the animal 
                                subfolders for a given time point. If no path is 
                                given the user is prompted to select the folder.
                                
    animal_name (str): This specifis the name of an individual animal to be 
                        analyzed. This changes the folder being looped over to 
                        the animal folder. If not provided all folders are
                        processed.
                        
    Returns:
    
    im_channel_paths (list): List containing the paths to all the image channels
                            for each animal. The order is Macro, Vessel, NP, 
                            tissue_outline.
    
    stack_name_list (list): A list containing all the stack names (animal folder 
                            names)                    
    
    '''
        
    if top_level_folder_path == "":
        top_path = get_folder_path()  
    else:
        top_path = top_level_folder_path
        
    if not animal_name == "":
        #signal image
        top_path = f'{top_level_folder_path}\\{animal_name}'
        #remember to exclude folders starting with '_'       
    else:
        #generates list of all the folders in the "time folder" also known as 
        #im stack name folders

        #Loop of folder
        stack_name_list = []
        im_channel_paths = []
        im_stack_name_subfolder_list = next(os.walk(top_path))[1] 
        
        for folder in im_stack_name_subfolder_list:
            if folder.startswith('_'):
                continue
            else:
                folder_name= folder
                stack_name_list.append(folder_name)
                
                path = f'{top_path}\\{folder_name}'
                image_channels = identify_image_channel_paths(path, folder_name)
                im_channel_paths.append(image_channels)
 
    return(im_channel_paths, stack_name_list)    