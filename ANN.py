# !/usr/bin/env python
# -*-coding:utf-8 -*-

# File       : ANN.py
# Description：build ANN model


import os
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler 
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import accuracy_score, roc_auc_score, f1_score
from sklearn.metrics import classification_report,confusion_matrix
import data_integration
import numpy as np
from data_integration import qm_descriptors

def scan_parameters(x,y):
    
    x = np.delete(x, len(qm_descriptors), axis=1)
    
    hidden_layer_sizes = [(50, 50), (100,100), (200,200), (400,400)]    
    solver  = ['lbfgs', 'sgd', 'adam']
    alpha = [0.001,0.0001, 0.00001]
    
    params_dict = {'hidden_layer_sizes': hidden_layer_sizes, 'solver' : solver, 'alpha': alpha, 'max_iter': [400]}
    
    scaler = StandardScaler(copy=True, with_mean=True, with_std=True)
    scaler.fit(x) 
    x = scaler.transform(x) 
    
    grid_search = GridSearchCV(estimator = MLPClassifier(), param_grid = params_dict, cv = 3,  n_jobs = 12)
    grid_search.fit(x, y)
    print(grid_search.best_params_)
    
    return grid_search.best_params_


def train_NN(X_train, y_train, X_valid, y_valid, param_list = None):
    
    x_train = np.delete(X_train, len(qm_descriptors), axis=1)
    x_valid = np.delete(X_valid, len(qm_descriptors), axis=1)
    
    scaler = StandardScaler(copy=True, with_mean=True, with_std=True)
    scaler.fit(x_train)
    
    x_train = scaler.transform(x_train)  
    
    if param_list == None:
        nn_classifier = MLPClassifier(max_iter = 400)
        
    else:
        nn_classifier = MLPClassifier(solver=param_list['solver'], hidden_layer_sizes = param_list['hidden_layer_sizes'], alpha = param_list['alpha'])
 
    nn_classifier.fit(x_train, y_train)
    
    
    x_valid = scaler.transform(x_valid)  
    predictions = nn_classifier.predict(x_valid)
    #probabilities = clf.predict_proba(x_valid) 
    
        
    #print(confusion_matrix(y_valid,predictions))
    print(classification_report(y_valid,predictions)) 
    accuracy = (100*accuracy_score(y_valid, predictions))    
    
    return accuracy

def build_NN_classifier(filename, option, model_name = None):  
     
    
    # LOAD DATA
    descriptors = qm_descriptors
    X, Y = data_preprocess.load_data(filename, descriptors)
    
    # IF DOWNSAMPLING:
    #print('>> Down sampling.')
    #smaller_x, smaller_y = data_preprocess.do_down_sampling(X,Y)
    
    
    if option == 'default':
        
        print('Start')
        print('*-----------------------------*')
        print('1.Train default parameters.')
        
        accuracies_default = []
        for i in range(10):
            x_train, x_valid, y_train, y_valid = data_preprocess.split_data (X, Y, partition=0.20)
            accuracies_default.append(train_NN(x_train, y_train, x_valid, y_valid))
        
        print('three times accuracy: %.2f' % np.mean(accuracies_default))
        
    elif option == 'train':
        
        print('*************************')
        print('Find best parameters.')
        
        params = []
        accuracies = []

        for i in range(10):
            x_train, x_valid, y_train, y_valid = data_preprocess.split_data (X, Y, partition=0.20)
            best_parameters = scan_parameters(x_train, y_train)
            params.append(best_parameters)
            accuracy = train_NN(x_train, y_train, x_valid, y_valid, best_parameters)
            accuracies.append(accuracy)
            
        print('**********************************')
        print('End.')
        print('**********************************')
        
        
        for i in range(len(accuracies)):
            print('Run ' + str (i+1)+ ' ', params[i], ' : ', accuracies[i])

        

    elif option == 'test':
        
        print('Start')
        print('********************************')
        
        hidden_layer_sizes = (100, 100)
        solver  = 'adam'
        alpha =  0.001
        
        params_dict = {'hidden_layer_sizes': hidden_layer_sizes, 'solver' : solver, 'alpha': alpha,  'max_iter': [400]}
    
 
        print(params_dict)
        
        acc_list = []
        for i in range(10):
            x_train, x_valid, y_train, y_valid = data_preprocess.split_data (X, Y, partition=0.20)
            acc_list.append(train_NN(x_train, y_train, x_valid, y_valid,params_dict))
        
        print('Result.')
        print('**********************************')
        print('Accuracy: %.2f' % np.mean(acc_list))
        

 

        


