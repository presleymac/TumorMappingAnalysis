import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import re

from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn import metrics
from sklearn.model_selection import cross_val_score

from imblearn.over_sampling import SMOTE   
from imblearn.over_sampling import SMOTENC   

from sklearn import feature_selection
from  sklearn.feature_selection import RFE
import statsmodels.api as sm
from sklearn.metrics import roc_curve
import neighbourhood_analysis_functions as nf

import itertools
# %matplotlib inline


def correlation(dataframe, threshold): 
    '''correlation
    
    Identifies sets of features that are correlated more than than threshod value
    and adds one of them to a set in an interative process. 
    
    Arguments: 
    
    dataframe(pd.DataFrame): DataFrame containing all features
    
    threshold(int):  Correlation threshold
    
    Returns:
    
    correlated_columns(set):  Set containing correlated columns
    '''
    correlated_columns = set()
    corr_matrix = dataframe.corr()

    for i in range(len(corr_matrix.columns)):
        for j in range(i):
            if abs(corr_matrix.iloc[i,j]) > threshold:
                colname = corr_matrix.columns[i]
                correlated_columns.add(colname)

    return correlated_columns


def load_single_radius_df(top_folder, radius):
    ''' load_single_radius_df
    load radius feature dataframe
    
    Arguments" 
    
    top_folder(str): path to top level folder
    radius(int): Nanoparticle region radius to process
    
    Returns: 
    radius_dataframe(pd.DataFrame): Dataframe coming single radius feature data
    
    '''
    
    data_folder = f'{top_folder}/Neighbourhood_Analysis_2023/Vessel_Analysis_Summary'
    file_name =  f'{data_folder}/{radius}um_radius_full_folder_vessel_data_um_2023.csv'
    radius_dataframe = pd.read_csv(file_name, index_col=[0])
    return(radius_dataframe)

def format_columns(data, drop_columns=['np_positive_region',
                                    'np_intensity',
                                    'ves_label',
                                    'total_macro_classified_pix',
                                    'num_hotspots_per_ves_seg', 
                                    'macro_count', 
                                    'vessel_count', 
                                    'total_ves_classified_pix', 
                                    'pre_classification_macro_pix',
                                    'pre_classification_ves_pix',
                                    'total_tissue_pixels']):
    
    '''removes non-numeric columns and drops the columns we don't want to work 
    with. Formats data for taining and texting'''
    
    #Find column types and drop any that are not numeric
    data.dtypes.unique()
    num = ['int64', 'float64']
    num_vars = list(data.select_dtypes(include=num))
    data = data[num_vars]
    data_copy = data.copy()
    
    X_data = data_copy.drop(columns=drop_columns, axis=1)
    y_data = data['np_positive_region']
    
    return(X_data, y_data)

def train_test_split_data(X,y, original_df, test_size=0.3, random_state=123, stratify=True):
    '''train_test_split_data
    
    Splits data into training and testing data stratified based on original 
    data category proportions
    
    Arguments:
    
    X(pd.DataFrame): features to be spit into training and testing data
    y(pd.DataFrame):  data categories
    original_df(pd.DataFrame): Original dataframe X and y are issolated from 
    test_size(int):  Ratio of testing to training data desired
    random_state(int):  Random seed default value 123. Any int can be used
    stratify(logical): default = True. Stratifys testing data such that proportion
                        of classes matches the original dataset
    
    Returns: 
    X_train(pd.DataFrame): X Training data
    y_train(pd.DataFrame): y training data
    X_test(pd.DataFrame): X testing data
    y_test(pd.DataFrame): y testing data
    
    '''
    
    starting_data = pd.concat([X,y], axis=1)
    if stratify == True:
        train, test = train_test_split(starting_data, 
                                    test_size = test_size, 
                                    random_state=random_state, 
                                    stratify=starting_data.np_positive_region)
    else:
        train, test = train_test_split(starting_data, 
                                    test_size = test_size, 
                                    random_state=random_state)
                
    y_train = train['np_positive_region']
    X_train = train.drop(['np_positive_region'], axis=1)
    y_test = test['np_positive_region']
    X_test = test.drop(['np_positive_region'], axis=1)
    return(X_train, y_train, X_test, y_test)

def display_correlation_matrix(X_train):
    '''display correlation matrix
    Displays heatmap of correlations between features within the supplied
    dataframe
    
    X_train(pd.DataFrame): Dataframe to show heatmap for
    
    '''
    
    corrmatrix= X_train.corr()
    sns.heatmap(corrmatrix)

def normalize_and_scale_data_smote_nc(X_train, y_train, X_test, y_test, top_folder_path, radius):
    '''normalize_and_scale_data_smote_nc
    
    Use SMTOE-NC to balance and then save datasets.
    
    Arguments:
    
    X_train(pd.DataFrame): Dataframe containing the training dataset
    y_train(pd.DataFrame): Dataframe contianing the training dataset classifications 
    X_test(pd.DataFrame): Dataframe containing testing dataset
    y_test(pd.DataFrame): Testing dataset corresponding classifications
    top_folder_path(str): Path to the folder top folder (folder level 1 in read me)
    radius(int): Integer value corresponding to the nanopartcile region radius size 
                in um
    
    Returns:
    normalized_X_train(pd.DataFrame): Training dataset after SMOTE-NC and normalization
    y_resampled(pd.DataFrame): Training y data after SMOTE-NC
    scaled_X_test(pd.DataFrame): Testing dataset after normalization
    
    
    Note: smote_nc is better for this analysis because have categorical variables 
    like time and np_size but smote_nc has a negligible/ no impact on the resultant 
    logit and takes significantly longer to run
    '''
    
    scaler = StandardScaler()
    smote_nc = SMOTENC(categorical_features=[11,12], random_state=0)
    X_resampled, y_resampled = smote_nc.fit_resample(X_train, y_train)
    normalized_X_train= scaler.fit_transform(X_resampled)
    
    scaled_X_test = pd.DataFrame(scaler.transform(X_test), columns=X_test.columns)
    X_norm_df = pd.DataFrame(normalized_X_train, columns=X_train.columns)
    
    save_folder = f'{top_folder_path}/Neighbourhood_Analysis_2023/SMOTE_NC_DATA'
    save_name_xtrain_org = (f'{save_folder}/X_training_data_org_{radius}um_radius_full_tumour.csv')
    save_name_ytrain_org = (f'{save_folder}/y_training_org_{radius}um_radius_full_tumour.csv')
    save_name_xtest_org = (f'{save_folder}/X_testing_org_{radius}um_radius_full_tumour.csv')
    save_name_ytest_org = (f'{save_folder}/y_testing_org_{radius}um_radius_full_tumour.csv')
    X_train.to_csv(save_name_xtrain_org)
    y_train.to_csv(save_name_ytrain_org)
    X_test.to_csv(save_name_xtest_org)
    y_test.to_csv(save_name_ytest_org)
    
    save_name_xnorm = (f'{save_folder}/X_training_normalized_smote_nc_{radius}um_radius_full_tumour.csv')
    save_name_xtrain = (f'{save_folder}/X_training_resampled_data_smote_nc_sampled_{radius}um_radius_full_tumour.csv')

    save_name_ytrain = (f'{save_folder}/y_training_data_smote_nc_sampled_{radius}um_radius_full_tumour.csv')
    save_name_xtest = (f'{save_folder}/X_testing_data_scaled_smote_nc_sampled_{radius}um_radius_full_tumour.csv')
    save_name_ytest = (f'{save_folder}/y_testing_data_smote_nc_sampled_{radius}um_radius_full_tumour.csv')
    
    
    X_norm_df.to_csv(save_name_xnorm)
    X_resampled.to_csv(save_name_xtrain)
    # X_resampled_df = pd.DataFrame(X_resamed)
    y_resampled.to_csv(save_name_ytrain)
    scaled_X_test.to_csv(save_name_xtest)
    y_test.to_csv(save_name_ytest)
    

    return( normalized_X_train, y_resampled, scaled_X_test)

def normalize_data(X_resampled, y_resampled, X_test):
    '''normalize_data
    
    This function defines a standard scaler for the model training data 
    (after any SMOTE-NC to balance the data) and then applies that scaler 
    to the testing data
    
    
    Arguments:
    X_resampled(pd.DataFrame): X_training data after any resampling to 
                                balance the data classes
    y_resampled(pd.DataFrame): y_training data after resampling to balance
                                classes
    X_test(pd.DataFrame):  Testing dataset to apply scaler to.
    
    Returns:
    
    normalized_X_train(pd.DataFrame): Training dataframe after scaler is applied
    y_resampled(pd.DataFrame): Training dataset class labels (unaltered)
    scaled_X_test(pd.DataFrame):  Testing dataset after scaler is applied
    '''
    
    scaler = StandardScaler()
    normalized_X_train= scaler.fit_transform(X_resampled)
    scaled_X_test = pd.DataFrame(scaler.transform(X_test), columns=X_test.columns)
    
    return(normalized_X_train, y_resampled, scaled_X_test )

def normalize_and_scale_data_smote(X_train, y_train, X_test, y_test):
    '''normalize_and_scale_data_smote
    
    Runs regular SMOTE sampling and then applies scaler
    
    Arguments:
    X_train(pd.DataFrame): X_training data
    y_train(pd.DataFrame): y_training data 
    X_test(pd.DataFrame):  Testing dataset to apply scaler to.
    y_test(pd.DataFrame): y testing dataset
    
    Returns:
    
    normalized_X_train(pd.DataFrame): SMOTE Training dataframe after scaler is applied
    y_resampled(pd.DataFrame): Training dataset after SMOTE
    scaled_X_test(pd.DataFrame):  Testing dataset after scaler is applied
    '''
    
    scaler = StandardScaler()
  
    smote = SMOTE(random_state=0)
    X_resampled, y_resampled = smote.fit_resample(X_train, y_train)
    normalized_X_train= scaler.fit_transform(X_resampled)
    scaled_X_test = pd.DataFrame(scaler.transform(X_test), columns=X_train.columns)
    # normalized_X_train = pd.DataFrame(normalized_X_train, columns=X_train.columns)
    return(normalized_X_train, y_resampled, scaled_X_test)

def prepare_training_data_for_logit(resampled_x_data, resampled_y_data, X_train):
    '''prepare_training_data_for_logit
    
    Prepare training data for logisitic regression ananlysis by removing
    colinear features with a pearson correlation greater than 0.7 and 
    add a constant
    
    Arguments: 
    resampled_x_data(pd.DataFrame): SMOTE-NC and normalized X training dataframe
    resampled_y_data(pd.DataFrame): SMOTE-NC and y training dataframe
    X_train(pd.DataFrame): Orginal X training dataset
    
    Returns: 
    X_train_intercept(pd.DataFrame): Finalized training dataset after dropping
                                    colinear features and adding 
                                    logit intercept (constant)
    corr_feature(set): Set containing correlated features withing training dataset
    
    
    '''
    
    
    X_resampled_df = pd.DataFrame(resampled_x_data, 
                                  columns= X_train.columns, 
                                  index=resampled_y_data.index)
    
    corr_feature = correlation(X_resampled_df, 0.7)
    X_resampled_df.drop(corr_feature, axis=1, inplace=True)
    X_train_intercept = sm.add_constant(X_resampled_df)
    return(X_train_intercept, corr_feature)


def run_logistic_regression_statsmodels(X_training_data_w_intercept, y_training_data):
    '''run_logistic_regression_statsmodels
    
    Train and run logistic regression using the stats models package
    
    Arguments:
    
    X_training_data_w_intercept(pd.DataFrame): Fully pre-processed logitic regression data
    y_training_data(pd.DataFrame): Corresponding classificaion labels for X data
    
    Returns:
    logistic_reg: Fitted statsmodel logisitic regression model
    training_cm(): Training data confusion matrix
    training_accuracy(int): Prediction accuracy on testing dataset
    
    '''
    
    logistic_reg =  sm.Logit(y_training_data, X_training_data_w_intercept).fit()
    print(logistic_reg.summary())
    
    training_predictions = np.round(logistic_reg.predict(X_training_data_w_intercept))
    training_cm = confusion_matrix(y_training_data, training_predictions)                         
    print ("Confusion Matrix : \n", training_cm)     
    training_accuracy = accuracy_score(y_training_data, training_predictions)

    return(logistic_reg, training_cm, training_accuracy)

def prepare_testing_data (scaled_X_test, training_correlated_features):
    '''prepare_testing_data 
       
       Prepare testing data for logit by dropping correlated features 
       dropped in the training dataset and adding an intercept
        
        Arguments:
        
        scaled_X_test(pd.DataFrame):  X testing data after scaler is applied
        training_correlated_features(list): List of correlated features from 
                                            training dataset
                                            
        Returns:
        pre_processed_X_test(pd.DataFrame): x_test data after preprocessing is preformed
    
    '''
    
    
    
    
    scaled_X_test.drop(training_correlated_features, axis=1, inplace=True)
    scaled_X_test = sm.add_constant(scaled_X_test)
    pre_processed_X_test =  scaled_X_test
    return(pre_processed_X_test)

def make_predictions(sm_logisitc_regression_model, pre_processed_X_test, y_test):
    '''make_predictions
    
    Make predictions with new data using the fitted statsmodels logistic regression model
    
    
    Arguments:
    sm_logistic_regression_model: fitted statsmodels logistic regression model
    
    pre_processsed_X_test(pd.Dataframe): Testing dataset with the same inputs as the training data
                                        Data must be transformed with the scaler fit to the
                                        training data.
    
    y_test(pd.DataFrame): Corresponding classification for features in the testing dataset. Used to 
            calculate testing accuracy
            
            
    Returns:
    testing_data_cm: Testing data Confusion matrix
    testing_accuracy(int): testing data prediction accuracy
    y_predictions(np.ndarray): model predictions for testing data
    
    '''
    
    
    y_predictions = np.round(sm_logisitc_regression_model.predict(pre_processed_X_test))
    testing_data_cm = confusion_matrix(y_test, y_predictions) 
    testing_accuracy = accuracy_score(y_test, y_predictions)
    return(testing_data_cm, testing_accuracy, y_predictions)

def scikit_learn_logit_modeling(X_training_data_w_intercept, 
                                y_training_data, 
                                pre_processed_X_test,
                                y_test):
    '''fit_int set as false because we have added a constant column without this
    we would be adding a new intercept column'''
    sk_logreg_intercept = LogisticRegression(fit_intercept=False)
    sk_logreg_intercept.fit(X_training_data_w_intercept, y_training_data)
    sk_y_pred_intercept = sk_logreg_intercept.predict(pre_processed_X_test)
    sk_y_pred_prob_intercept =  sk_logreg_intercept.predict_proba(pre_processed_X_test)[:,1]
    
    
    auc = metrics.roc_auc_score(y_test, sk_y_pred_prob_intercept)
    print('AUC = ',auc)
    #false_positives, true_positives
    fpr, tpr, thresholds=  roc_curve(y_test, sk_y_pred_prob_intercept)
    plt.plot(fpr, tpr, label = 'Logisitic Regression')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Logisitic Regression ROC Curve')
    plt.show()
    
    sk_logit_cm = confusion_matrix(y_test, sk_y_pred_intercept) 
  
    sk_acc = accuracy_score(y_test, sk_y_pred_intercept)
    
    return(sk_logreg_intercept, sk_y_pred_intercept, sk_y_pred_prob_intercept, fpr, tpr, thresholds, auc, sk_acc)


def cross_validate_sklearn_logistic_reg(sklearn_fitted_logistic_reg_model,
                                        X_training_data_with_interecept,
                                        y_training_data,
                                        cross_validate_number =5):
    
    accuracry_cross_val_scores = cross_val_score(sklearn_fitted_logistic_reg_model, 
                                                X_training_data_with_interecept, 
                                                y_training_data, 
                                                cv=cross_validate_number)
    
    f1_cross_val_scores = cross_val_score(sklearn_fitted_logistic_reg_model, X_training_data_with_interecept, y_training_data, cv=cross_validate_number, scoring='f1')
    roc_auc_cross_val_scores = cross_val_score(sklearn_fitted_logistic_reg_model, 
                                                X_training_data_with_interecept, 
                                                y_training_data, 
                                                cv=cross_validate_number, 
                                                scoring='roc_auc')
                        
    return(accuracry_cross_val_scores, f1_cross_val_scores, roc_auc_cross_val_scores)


def powerset(iterable):
    s = list(iterable)
    return (itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(1, len(s)+1)))

def statsmodels_logistic_regression_optimization(X_training_data_w_intercept, y_training_data):
    count = 0
    features = X_training_data_w_intercept.columns
    
    aic_list = []
    bic_list = []

    llf_list = []
    llr_list= []
    # llr_pval_list = []
    prsquared_list = []
    feature_list = [] 
    for feature_subset in powerset(features):
        #   print(subset)
        subset_data = (X_training_data_w_intercept[list(feature_subset)])
        
        count += 1
        sm_model =  sm.Logit(y_training_data, subset_data).fit()
        aic = sm_model.aic
        bic = sm_model.bic
        llr = sm_model.llr
        llf = sm_model.llf
        prsq = sm_model.prsquared  
  
        aic_list.append(aic)
        bic_list.append(bic)
        feature_list.append(feature_subset)
        llf_list.append(llf)
        llr_list.append(llr)

        prsquared_list.append(prsq)
    
    optimization_df= pd.DataFrame(columns= ['Features',
                                            'AIC', 
                                            'BIC', 
                                            'Likelihood Ratio Chi-Squared Statistic', 
                                            'Log Likelihood Fit', 
                                            'Pseudo R-squared'] )  
    
    optimization_df['Features'] = feature_list
    optimization_df['AIC'] = aic_list
    optimization_df['BIC'] = bic_list
    optimization_df['Likelihood Ratio Chi-Squared Statistic'] = llr_list
    optimization_df['Log Likelihood Fit'] = llf_list
    optimization_df['Pseudo R-squared'] = prsquared_list
    
    
    best_aic = (np.where(aic_list ==min(aic_list) ))
    best_bic = (np.where(bic_list == min(bic_list)))
    best_llf = (np.where(llf_list==max(llf_list )))
    best_lr = (np.where(llr_list==max(llr_list)))
    best_prsq = (np.where(prsquared_list==(max(prsquared_list))))
    
    best_df = pd.DataFrame(columns=['Features Best Pseudo R-sq',
                            'AIC', 
                            'BIC', 
                            'Likelihood Ratio Chi-Squared Statistic', 
                            'Log Likelihood Fit', 
                            'Pseudo R-squared'])
    
    best_df['Features Best Pseudo R-sq'] = [feature_list[best_prsq[0][0]]]
    best_df['AIC'] = best_aic
    best_df['BIC'] = best_bic
    best_df['Likelihood Ratio Chi-Squared Statistic'] = best_lr
    best_df['Log Likelihood Fit'] = best_llf
    best_df['Pseudo R-squared'] = best_prsq

    return(optimization_df, best_df)
                                                
def save_optimaztion_df(top_path, radius, optimization_df, best_df):
    save_folder = f'{top_path}/Neighbourhood_Analysis_2023/Summary_Dataframe'
    save_name = (f'{save_folder}/Logistic_Regression_optimization_{radius}um_radius_full_tumour.csv')
    save_optimal = (f'{save_folder}/Logistic_Regression_optimization_{radius}um_radius_optimal_featuers_full_tumour.csv')
    optimization_df.to_csv(save_name)
    best_df.to_csv(save_optimal)
    
def run_logistic_reg_smote_NC(top_folder_path, radius):   
    data = load_single_radius_df(top_folder_path, radius)
    data.pivot_table(index='np_positive_region', aggfunc='size').plot(kind='bar')
    X, y = format_columns(data)
    X_train, y_train, X_test, y_test = train_test_split_data(X, y, data)
    display_correlation_matrix(X_train)
    
    X_resampled, y_resampled, scaled_X_test =normalize_and_scale_data_smote_nc(X_train, y_train, X_test, y_test)
    correlated_features = correlation(X_resampled, 0.7)
    
    save_folder = f'{top_path}/Neighbourhood_Analysis_2023/SMOTE_NC_DATA'
    
    save_name_xtrain_org = (f'{save_folder}/X_training_data_org{radius}um_radius_full_tumour.csv')
    save_name_ytrain_org = (f'{save_folder}/y_training_org{radius}um_radius_full_tumour.csv')
    save_name_xtest_org = (f'{save_folder}/X_testing_org{radius}um_radius_full_tumour.csv')
    save_name_ytest_org = (f'{save_folder}/y_testing_org{radius}um_radius_full_tumour.csv')
    
    X_train.to_csv(save_name_xtrain_org)
    y_train.to_csv(save_name_ytrain_org)
    X_test.to_csv(save_name_xtest_org)
    y_test.to_csv(save_name_ytest_org)
    
    save_name_xtrain = (f'{save_folder}/X_training_data_smote_nc_sampled{radius}um_radius_full_tumour.csv')
    save_name_ytrain = (f'{save_folder}/y_training_data_smote_nc_sampled{radius}um_radius_full_tumour.csv')
    save_name_xtest = (f'{save_folder}/X_testing_data_smote_nc_sampled{radius}um_radius_full_tumour.csv')
    save_name_ytest = (f'{save_folder}/y_testing_data_smote_nc_sampled{radius}um_radius_full_tumour.csv')
    
    
    X_train_w_intercept, corr_feature =prepare_training_data_for_logit(X_resampled, y_resampled, X_train)
    X_train_w_intercept=X_train_w_intercept.drop(columns=['num_hotspots_per_ves_seg'], axis=1)
    scaled_X_test = scaled_X_test.drop(columns=['num_hotspots_per_ves_seg'], axis=1)
    
    
    
    logistic_reg, training_cm, training_accuracy = run_logistic_regression_statsmodels(X_train_w_intercept, y_train)
    pre_processed_X_test = prepare_testing_data (scaled_X_test, correlated_features)
    
    testing_data_cm, testing_accuracy =  make_predictions(logistic_reg, pre_processed_X_test, y_test)
    
    optimization_df, best_df = statsmodels_logistic_regression_optimization(X_training_data_w_intercept,  y_training_data)
    
    
    return(X_train_w_intercept, y_resampled, scaled_X_test, pre_processed_X_test,optimization_df,best_df)    
    

    
if __name__ == '__main__':
    top_folder_path = nf.get_folder_path()
    radius_list = [25]
    for radius in radius_list:
        print(radius)
        data = load_single_radius_df(top_folder_path, radius)
        data.pivot_table(index='np_positive_region', aggfunc='size').plot(kind='bar')
        
        
        X, y = format_columns(data)
        X_train, y_train, X_test, y_test = train_test_split_data(X, y, data)
        display_correlation_matrix(X_train)
    
        
        # X_resampled, y_resampled, scaled_X_test =normalize_and_scale_data_smote(X_train, y_train, X_test, y_test)
        # correlated_features = correlation(X_resampled, 0.7)
        
        X_norm, y_resampled, scaled_X_test =normalize_and_scale_data_smote_nc(X_train, y_train, X_test, y_test, top_folder_path, radius)

        # correlated_features = correlation(X_resampled, 0.7)


        
        X_train_w_intercept, corr_feature=prepare_training_data_for_logit(X_norm, y_resampled, X_train)
        
        
        # X_train_w_intercept=X_train_w_intercept.drop(columns=['num_hotspots_per_ves_seg'], axis=1)
        # scaled_X_test = scaled_X_test.drop(columns=['num_hotspots_per_ves_seg'], axis=1)
        # X_train_w_intercept.to_csv(save_name_xtrain)
        # X_train_w_intercept=prepare_training_data_for_logit(X_resampled, y_resampled, X_train)
        
        logistic_reg, training_cm, training_accuracy = run_logistic_regression_statsmodels(X_train_w_intercept, y_resampled)
        # correlated_features = correlation(X_resampled, 0.7)
        pre_processed_X_test = prepare_testing_data (scaled_X_test, corr_feature)
        
        testing_data_cm, testing_accuracy, y_pred_logit =  make_predictions(logistic_reg, pre_processed_X_test, y_test)
        
        optimization_df, best_df = statsmodels_logistic_regression_optimization(X_train_w_intercept,  y_resampled)
        
        
        save_optimaztion_df(top_folder_path, radius, optimization_df, best_df)

        best_features = best_df['Features Best Pseudo R-sq'][0]
        best_model_X_features= X_train_w_intercept[list(best_features)]
        
        
        best_logistic_reg, best_training_cm, best_training_accuracy = run_logistic_regression_statsmodels(best_model_X_features, y_resampled)
        best_logit_feature_results = pd.DataFrame(best_logistic_reg.summary2().tables[1])
        best_logit_summary_results = pd.DataFrame(best_logistic_reg.summary2().tables[0])
        testing_data_cm_best, testing_accuracy_best, y_pred_best =  make_predictions(best_logistic_reg, pre_processed_X_test, y_test)

        
        
        best_logit_feature_results.to_csv(f'{top_folder_path}/Neighbourhood_Analysis_2023/'
                                        f'Logistic_Regression_Results/'
                                        f'Optimized_Model_SMOTE_NC_{radius}um_radius'
                                        f'_feature_results.csv')

        
        best_logit_summary_results.to_csv(f'{top_folder_path}/Neighbourhood_Analysis_2023/'
                                        f'Logistic_Regression_Results/'
                                        f'Optimized_Model_SMOTE_NC_{radius}um_radius'
                                        f'_summary_results.csv')

        sk_logreg, sk_y_pred, sk_y_pred_prob, fpr, tpr, thresholds, auc, sk_acc = scikit_learn_logit_modeling(X_train_w_intercept,
                                                                                                            y_resampled, 
                                                                                                            pre_processed_X_test,
                                                                                                            y_test)
        
        print(classification_report(y_test, sk_y_pred))
        
        accuracry_cross_val_scores, f1_cross_val_scores, roc_auc_cross_val_scores= cross_validate_sklearn_logistic_reg(sk_logreg,
                                                                                                                    X_train_w_intercept,
                                                                                                                    y_resampled,
                                                                                                                    cross_validate_number =5)
        # X_train, y_train, X_test, y_test = train_test_split_data(X, y, data)
        # display_correlation_matrix(X_train)
        # X_resampled, y_resampled, scaler = normalize_data_smote_nc(X_train, y_train)
        
        # X_resampled_df = pd.DataFrame(X_resampled, 
        #                               columns= X_train.columns, 
        #                               index=y_resampled.index)
        
        # corr_feature = correlation(X_resampled, 0.7)
        # X_resampled_df.drop(corr_feature, axis=1, inplace=True)
        # X_train_intercept = sm.add_constant(X_resampled_df)
        

        report = classification_report(y_test, sk_y_pred, output_dict=True)
        classification_df = pd.DataFrame(report).transpose()
        classification_df.to_csv(f'{top_folder_path}/Neighbourhood_Analysis_2023/'
                                        f'Logistic_Regression_Results/{radius}um_radius'
                                        f'_classification_report.csv')
        training_cm_df = pd.DataFrame(training_cm)
        training_cm_df.to_csv(f'{top_folder_path}/Neighbourhood_Analysis_2023/'
                                        f'Logistic_Regression_Results/{radius}um_radius'
                                        f'_training_confusion_matrix.csv')
        testing_data_cm_df = pd.DataFrame(testing_data_cm)
        testing_data_cm_df.to_csv(f'{top_folder_path}/Neighbourhood_Analysis_2023/'
                                        f'Logistic_Regression_Results/{radius}um_radius'
                                        f'_testing_confusion_matrix.csv')