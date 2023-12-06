import os
import numpy as np 
import pandas as pd 
from skimage import io 
import time
import neighbourhood_analysis_functions as nf 





def index_of_dispersion(dataframe, dataframe_column):
    '''Calculates the index of dispersion for a given column in a dataframe'''
    data = dataframe[dataframe_column]
    varience = np.var(data)
    mean = np.mean(data)
    IOD = varience/mean
    return(IOD)


def coef_variation (dataframe, dataframe_column):
    '''Calculates the coefficient of variation for a given column in a dataframe'''
    data = dataframe[dataframe_column]
    stdev = np.std(data)
    mean = np.mean(data)
    coef = stdev/mean
    return(coef)


def index_of_dispersion_image_folder_loop(top_folder, sample_number, radius_list):
    '''Calculates IOD for a given sample number for each image subfolder for all radius sizes'''
    
    
    #From the top level folder generate a list of all the image subfolder to analyze
    
    image_subfolder_list = nf.generate_list_subfolders(top_folder)
    msc_list = [image for image in image_subfolder_list if 'MSC' in image ]
    ut_list = [image for image in image_subfolder_list if 'UT' in image ]
    full_image_subfolder_list = msc_list + ut_list
    
    #Create paths and folders to store data
    neighbourhood_analysis_full_folder_path = f'{top_folder}/Neighbourhood_Analysis_2023'
    if not os.path.exists(neighbourhood_analysis_full_folder_path):
        os.mkdir(neighbourhood_analysis_full_folder_path)
    
    IOD_folder_path = f'{neighbourhood_analysis_full_folder_path}/IOD_Summary_DataFrames'
    if not os.path.exists(IOD_folder_path):
        os.mkdir(IOD_folder_path)
    

    time_str = top_folder.split('\\')[-1]
    time_pt = float((time_str).split('h')[0])
    tum_np_size = top_folder.split('\\')[-2]
    size_str = tum_np_size.split("GNP")[-1]
    np_size =  size_str.split('nm')[0]
    
    for radius in radius_list:  
        IOD_folder_level_save_path = f'{IOD_folder_path}/{tum_np_size}_{time_str}_IOD_um_calculations.csv'
  
        for image_path in full_image_subfolder_list:
            
            time_list = []
            animal_name_list = []
            image_name_list = []
            IOD_macro_list = []
            coef_of_var_macro_list = []
            IOD_ves_list = [] 
            coef_of_var_ves_list =[]
            IOD_np_list = []
            coef_of_var_np_list = [] 
            radius_size_list = [] 
            sample_number_list = []
            np_size_list = []
            
            IOD_df = pd.DataFrame()
            save_path = f'{image_path}/Neighbourhood_Analysis_2023/IOD_DataFrames'
            summary_df_path = f'{image_path}/Neighbourhood_Analysis_2023/Summary_Dataframes'
            if not os.path.exists(save_path):
                os.mkdir(save_path)

            stack_name =  nf.return_stack_name(image_path)   
            iod_sampling_data = (f'{summary_df_path}/{stack_name}_{sample_number}'
                                f'x_{radius}um_sphere_sampling_index_of_dispersion.csv')
            iod_unprocessed_data = pd.read_csv(iod_sampling_data)
            #Drops data less than the 10% tissue threshold
            drop_data = np.where(iod_unprocessed_data['pct_tissue']<10)[0]
            iod_unprocessed_data = iod_unprocessed_data.drop(drop_data)
            columns = ['pct_macro_per_tissue', 'pct_ves_per_tissue', 'pct_np_per_tissue' ]
            size = radius    
            iod_macro = index_of_dispersion(iod_unprocessed_data, columns[0])
            iod_ves= index_of_dispersion(iod_unprocessed_data, columns[1])
            iod_np = index_of_dispersion(iod_unprocessed_data, columns[2])
            cov_macro = coef_variation(iod_unprocessed_data, columns[0])
            cov_ves= coef_variation(iod_unprocessed_data, columns[1])
            cov_np = coef_variation(iod_unprocessed_data, columns[2])
            animal_name = stack_name.split('-')[0]
            
            animal_name_list.append(animal_name)
            image_name_list.append(stack_name)
            IOD_macro_list.append(iod_macro)
            coef_of_var_macro_list.append(cov_macro)
            IOD_ves_list.append(iod_ves) 
            coef_of_var_ves_list.append(cov_ves)
            IOD_np_list.append(iod_np)
            coef_of_var_np_list.append(cov_np)
            radius_size_list.append(radius)
            time_list.append(time_pt)
            sample_number_list.append(len(iod_unprocessed_data))
            np_size_list.append(np_size)
            single_image_iod_path  = f'{save_path}/{stack_name}_{sample_number}_samples_{radius}um_radius_IOD_um_calculations.csv'

            
            IOD_df['Animal Name'] = animal_name_list
            IOD_df['Image Name'] = image_name_list
            IOD_df['IOD Macrophages'] = IOD_macro_list
            IOD_df['IOD Vessel'] = IOD_ves_list
            IOD_df['IOD Nanoparticles'] = IOD_np_list
            IOD_df['Coefficient of Variation Macrophages'] = coef_of_var_macro_list
            IOD_df['Coefficient of Variation Vessels'] = coef_of_var_ves_list
            IOD_df['Coefficient of Variation Nanoparticles'] = coef_of_var_np_list 
            IOD_df['Radius (um)'] = radius_size_list
            IOD_df['Time Point'] = time_list
            IOD_df['Nanoparticle Size'] = np_size_list
            IOD_df['Number of Samples'] = sample_number_list
            IOD_df.to_csv(single_image_iod_path)


def index_of_dispersion_full_time_pt(top_folder, sample_number, radius_list):
    '''Calculates IOD for a given sample number for each time point for all radius sizes'''
        
    image_subfolder_list = nf.generate_list_subfolders(top_folder)
    msc_list = [image for image in image_subfolder_list if 'MSC' in image ]
    ut_list = [image for image in image_subfolder_list if 'UT' in image ]
    full_image_subfolder_list = msc_list + ut_list
    
    neighbourhood_analysis_full_folder_path = f'{top_folder}/Neighbourhood_Analysis_2023'
    
    IOD_folder_path = f'{neighbourhood_analysis_full_folder_path}/IOD_Summary_DataFrames'
    if not os.path.exists(IOD_folder_path):
        os.mkdir(IOD_folder_path)
    
    time_str = top_folder.split('\\')[-1]
    # time_pt = float((time_str).split('h')[0])
    tum_np_size = top_folder.split('\\')[-2]
    
    # radius_df = pd.DataFrame()
    for radius in radius_list:
        IOD_folder_level_save_path = f'{IOD_folder_path}/{tum_np_size}_{time_str}_{sample_number}_samples_{radius}um_radius_IOD_um_calculations.csv'
        radius_df = pd.DataFrame()
        for image_path in full_image_subfolder_list:
            stack_name =  nf.return_stack_name(image_path)   
            neighbourhood_analysis_single_path = f'{top_folder}/{stack_name}/Neighbourhood_Analysis_2023'
            IOD_single_path  = f'{neighbourhood_analysis_single_path}/IOD_DataFrames'


            single_image_iod_path  = f'{IOD_single_path}/{stack_name}_{sample_number}_samples_{radius}um_radius_IOD_um_calculations.csv'
            image_df = pd.read_csv(single_image_iod_path, index_col=[0])
            radius_df = pd.concat([radius_df, image_df])

        radius_df.to_csv(IOD_folder_level_save_path)


def iod_single_size_combo_df(top_folder, sample_number, radius_list):
    '''Calculates IOD for a given sample number for each nanoparticle size for
     all radius sizes'''
    
    image_subfolder_list = nf.generate_list_subfolders(top_folder)
    msc_list = [image for image in image_subfolder_list if 'MSC' in image ]
    ut_list = [image for image in image_subfolder_list if 'UT' in image ]
    full_image_subfolder_list = msc_list + ut_list
    
    neighbourhood_analysis_full_folder_path = f'{top_folder}/Neighbourhood_Analysis_2023'
    
    IOD_folder_path = f'{neighbourhood_analysis_full_folder_path}/IOD_Summary_DataFrames'
    
    
    if not os.path.exists(IOD_folder_path):
        os.mkdir(IOD_folder_path)
    
    time_str = top_folder.split('\\')[-1]
    # time_pt = float((time_str).split('h')[0])
    tum_np_size = top_folder.split('\\')[-2]
    
    for image_path in full_image_subfolder_list:
        all_size_df = pd.DataFrame()
        stack_name =  nf.return_stack_name(image_path)   
        neighbourhood_analysis_single_path = f'{top_folder}/{stack_name}/Neighbourhood_Analysis_2023'
        IOD_single_path  = f'{neighbourhood_analysis_single_path}/IOD_DataFrames'
        IOD_save_path = f'{IOD_single_path}/{stack_name}_{sample_number}_samples_all_radius_size_IOD_um_units_data.csv'
        for radius in radius_list: 

            single_image_iod_path  = f'{IOD_single_path}/{stack_name}_{sample_number}_samples_{radius}um_radius_IOD_um_calculations.csv'
            image_df = pd.read_csv(single_image_iod_path, index_col=[0])
            all_size_df = pd.concat([all_size_df, image_df])

        all_size_df.to_csv(IOD_save_path)


def iod_all_sizes_all_images_folder_level_concat(top_folder, sample_number, radius_list):
    '''Calculates IOD for a given sample number for each all nanoparticle sizes 
     and for all radius sizes
    '''
        
    image_subfolder_list = nf.generate_list_subfolders(top_folder)
    msc_list = [image for image in image_subfolder_list if 'MSC' in image ]
    ut_list = [image for image in image_subfolder_list if 'UT' in image ]
    full_image_subfolder_list = msc_list + ut_list
    
    neighbourhood_analysis_full_folder_path = f'{top_folder}/Neighbourhood_Analysis_2023'
    
    IOD_folder_path = f'{neighbourhood_analysis_full_folder_path}/IOD_Summary_DataFrames'
    time_str = top_folder.split('\\')[-1]
    # time_pt = float((time_str).split('h')[0])
    tum_np_size = top_folder.split('\\')[-2]
    
    size_combo_df = pd.DataFrame()
    for radius in radius_list:
        IOD_folder_level_save_path = f'{IOD_folder_path}/{tum_np_size}_{time_str}_{sample_number}_samples_{radius}um_radius_IOD_um_calculations.csv'
        combo_save_str = f'{IOD_folder_path}/{tum_np_size}_{time_str}_{sample_number}_samples_all_radius_sizes_um_units.csv'
        single_size_df = pd.read_csv(IOD_folder_level_save_path, index_col=[0])
        size_combo_df = pd.concat([size_combo_df, single_size_df])
        size_combo_df.to_csv(combo_save_str)
        

def iod_single_radius_all_np_sizes_all_time_pts(top_folder, sample_number, radius_list):
    '''Calculates IOD for a given sample number for radius size for
     all nanopartcle sizes and time points'''
    
    np_size_folder_list = nf.generate_list_subfolders(top_folder)
    np_size_folder_list = [x for x in np_size_folder_list if (x.split('\\')[-1]).startswith('U87')]
    full_tumour_data = f'{top_folder}/Neighbourhood_Analysis_2023'
    if not os.path.exists(full_tumour_data):
        os.mkdir(full_tumour_data)
    

    for radius in radius_list:
        radius_df = pd.DataFrame()
        single_np_size_df = pd.DataFrame()
        save_path_all_sizes = f'{full_tumour_data}/{radius}um_radius_{sample_number}_samples_all_particle_sizes_um_units.csv'

        for np_size_folder in np_size_folder_list:
            np_size = np_size_folder.split('\\')[-1]
            analysis_path = f'{np_size_folder}/Neighbourhood_Analysis_2023'
            file_path = f'{analysis_path}/{np_size}_{sample_number}_samples_{radius}um_radius_all_time_pts_IOD_um_calculations.csv'
            np_size_df = pd.read_csv(file_path, index_col=[0])
            single_np_size_df = pd.concat([single_np_size_df, np_size_df])

            
        radius_df = pd.concat([radius_df, single_np_size_df])
        radius_df.to_csv(save_path_all_sizes)   


def iod_single_radius_all_time_pts_df(np_size_folder, sample_number, radius_list):
    '''Calculates IOD for a given sample number for a nanoparticle size for
     a single radius size across all time points'''
    
    
    image_subfolder_list = nf.generate_list_subfolders(np_size_folder)

    
    neighbourhood_analysis_full_folder_path = f'{np_size_folder}/Neighbourhood_Analysis_2023'
    
    IOD_folder_path = f'{neighbourhood_analysis_full_folder_path}/IOD_Summary_DataFrames'
    
    time_path_list = [path for path in image_subfolder_list if path.endswith('h')]


    
    if not os.path.exists(neighbourhood_analysis_full_folder_path):
        os.mkdir(neighbourhood_analysis_full_folder_path)
   
    if not os.path.exists(IOD_folder_path):
        os.mkdir(IOD_folder_path)
    tum_np_size = np_size_folder.split('\\')[-1]    
    for radius in radius_list:
        single_radius_save_path = f'{neighbourhood_analysis_full_folder_path}/{tum_np_size}_{sample_number}_samples_{radius}um_radius_all_time_pts_IOD_um_calculations.csv'
        radius_df = pd.DataFrame()
        for time_folder_path in time_path_list:
            time_str = time_folder_path.split('\\')[-1]
            stack_list_path = nf.generate_list_subfolders(time_folder_path)
            msc_list = [image for image in stack_list_path if 'MSC' in image ]
            ut_list = [image for image in stack_list_path if 'UT' in image ]
            full_image_subfolder_list = msc_list + ut_list 
            time_df = pd.DataFrame()
            for image_path in full_image_subfolder_list:
                stack_name =  nf.return_stack_name(image_path)   
                neighbourhood_analysis_single_path = f'{np_size_folder}/{time_str}/{stack_name}/Neighbourhood_Analysis_2023'
                IOD_single_path  = f'{neighbourhood_analysis_single_path}/IOD_DataFrames'


                single_image_iod_path  = f'{IOD_single_path}/{stack_name}_{sample_number}_samples_{radius}um_radius_IOD_um_calculations.csv'
                image_df = pd.read_csv(single_image_iod_path, index_col=[0])
           
                time_df = pd.concat([time_df, image_df])
            radius_df = pd.concat([radius_df, time_df])
        radius_df.to_csv(single_radius_save_path)


def run_all_iod_calculations(top_folder, sample_number, radius_list ):
    ''' run_all_iod_calculations
    runs all index of dispersion calculations
    Note: assumes you select the top folder (folder level 1 in Read Me)
    
    Arguments
    top_folder(str): Path to the top level folder
    sample_number(int): Number of spherical sampels
    radius_list (list): List of radius sizes to calulate the IOD for
    
    '''

    
    
    print('Processed Folder: ', top_folder)
    
    #Find nanoparticle size folder list
    np_size_folder_list = nf.generate_list_subfolders(top_folder)
    np_size_folder_list = [x for x in np_size_folder_list if (x.split('\\')[-1]).startswith('U87')]
    
    time_folder_list = []
    
    for size_path in np_size_folder_list:
        #Generate list of nanoparticle sizes
        size_list = nf.generate_list_subfolders(size_path)
    
        time_folder_list.extend(size_list)

    #Generate time point list
    time_folder_list = [x for x in time_folder_list if (x.split('\\')[-1]).endswith('h')]
    
    
    for time_pt_folder in time_folder_list:
        #run various IOD calculations and saves dataframes
        index_of_dispersion_image_folder_loop(time_pt_folder, sample_number, radius_list)
        index_of_dispersion_full_time_pt(time_pt_folder, sample_number, radius_list)
        iod_single_size_combo_df(time_pt_folder, sample_number, radius_list)
        iod_all_sizes_all_images_folder_level_concat(time_pt_folder, sample_number, radius_list)
        
    np_size_list = nf.generate_list_subfolders(top_folder)
    
    for np_size in np_size_list:
        iod_single_radius_all_time_pts_df(np_size, sample_number, radius_list)
    
    iod_single_radius_all_np_sizes_all_time_pts(top_folder, sample_number, radius_list)
        
def iod_summary_data_by_radius(top_folder, sample_number, radius_list):
    '''iod_summary_data_by_radius:
    Finds basical summary values for IOD based on tumour feature and saves them 
    a csv file 
    
    Arguments:
    top_folder(str): Path to top level folder
    sample_number(int): number of samples
    radius_list(list): List of radius values
    
    
    
    '''
    
    path_to_data = f'{top_folder}/Neighbourhood_Analysis_2023'
    for radius in radius_list:
        file_name = f'{path_to_data}/{radius}um_radius_{sample_number}_samples_all_particle_sizes_um_units.csv'
        iod_df = pd.read_csv(file_name, index_col=[0])
        summary_name = f'{path_to_data}/{radius}_um_radius_{sample_number}_samples_summary_measurements_um_units.csv'
        columns = iod_df.columns
        
        index_names = columns.drop(['Animal Name', 'Image Name'])
        summary_stats_df = pd.DataFrame(index=index_names, columns = ['Mean', 
                                                                      'Median', 
                                                                      'Min', 
                                                                      'Max',
                                                                      'Standard Deviation',
                                                                      'Sample Number',
                                                                      'Above Median Count', 
                                                                      'Above Median %', 
                                                                      "Above Mean Count", 
                                                                      'Above Mean %', 
                                                                      'IOD Less Than 1.5',
                                                                      'IOD Less Than 1.5 %',
                                                                      'IOD Less Than 2',
                                                                      'IOD Less Than 2 %'])
        for ind in index_names:
            mean = np.mean(iod_df[ind])
            med = np.median(iod_df[ind])
            col_min = np.min(iod_df[ind])
            col_max = np.max(iod_df[ind])
            col_std = np.std(iod_df[ind])
            length = len(iod_df[ind])
            above_med = np.sum(iod_df[ind]>med)
            above_med_pct = above_med/length*100
            above_mean = np.sum(iod_df[ind]>mean)
            above_mean_pct = above_mean/length*100
            under_15 = sum(iod_df[ind]<1.5)
            under_15_pct = under_15/length*100
            under_2 = sum(iod_df[ind]<2)
            under_2_pct = under_2/length*100

            summary_stats_df['Mean'].loc[ind] = mean 
            summary_stats_df['Median'].loc[ind]= med 
            summary_stats_df['Min'].loc[ind]= col_min
            summary_stats_df['Max'].loc[ind]= col_max
            summary_stats_df['Standard Deviation'].loc[ind]= col_std
            summary_stats_df['Sample Number'].loc[ind]= length
            summary_stats_df['Above Median Count'].loc[ind] = above_med
            summary_stats_df['Above Median %'].loc[ind] = above_med_pct
            summary_stats_df["Above Mean Count"].loc[ind] =above_mean
            summary_stats_df['Above Mean %'].loc[ind] =above_mean_pct
            summary_stats_df['IOD Less Than 1.5'].loc[ind] = under_15
            summary_stats_df['IOD Less Than 1.5 %'].loc[ind] = under_15_pct
            summary_stats_df['IOD Less Than 2'].loc[ind] = under_2
            summary_stats_df['IOD Less Than 2 %'].loc[ind] = under_2_pct
        summary_stats_df.to_csv(summary_name)
            
            
if __name__ == '__main__':
    top_folder_path = nf.get_folder_path()
    sample_number =1000
    radius_list = [25]
    image_subfolder_list = nf.generate_list_subfolders(top_folder_path)
    run_all_iod_calculations(top_folder_path, sample_number, radius_list)
    iod_summary_data_by_radius(top_folder_path, sample_number, radius_list)
    iod_single_radius_all_time_pts_df(top_folder_path, sample_number, radius_list)

    
