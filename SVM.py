# !/usr/bin/env python
# -*-coding:utf-8 -*-

# File       : SVM.py
# Description：build SVM model

from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV

from sklearn.metrics import accuracy_score, roc_auc_score, f1_score
from sklearn.metrics import classification_report, confusion_matrix

import data_integration
import numpy as np

from data_integration import qm_descriptors

from sklearn.feature_selection import RFE



def scan_parameters(x, y):
    """
        svm types are :
        linear  []
        poly    [requires to specify the degree of polynomial]
        rbf     [gaussian] --> default
        sigmoid [only for bi-classifier]
    """
    x = np.delete(x, len(qm_descriptors), axis=1)
    svr = SVC()
    kernels = ['linear', 'rbf', 'sigmoid']
    Cs = [0.001, 0.01, 0.1, 1, 10, 100, 200]
    gammas = ['auto', 0.001, 0.01, 0.1, 1, 10]
    weights = [None]
    params_dict = {'kernel': kernels, 'C': Cs, 'class_weight': weights, 'gamma': gammas}
    scaler = StandardScaler(copy=True, with_mean=True, with_std=True)
    scaler.fit(x)
    x = scaler.transform(x)
    grid_search = GridSearchCV(estimator=svr, param_grid=params_dict, cv=3, n_jobs=12)
    grid_search.fit(x, y)
    print(grid_search.best_score_, grid_search.best_params_)

    return grid_search.best_params_


def train_SVM(X_train, y_train, X_valid, y_valid, param_list=None):
    x_train =np.delete(X_train, len(qm_descriptors), axis=1)
    x_valid = np.delete(X_valid, len(qm_descriptors), axis=1)
    scaler = StandardScaler(copy=True, with_mean=True, with_std=True)
    scaler.fit(x_train)
    x_train = scaler.transform(x_train)
    if param_list == None:
        my_SVM_classifier = SVC()
    else:
        my_SVM_classifier = SVC(kernel=param_list['kernel'], C=param_list['C'], class_weight=param_list['class_weight'],
                                gamma=param_list['gamma'])

    my_SVM_classifier.fit(x_train, y_train)
    x_valid = scaler.transform(x_valid)
    predictions = my_SVM_classifier.predict(x_valid)
    # probabilities = clf.predict_proba(x_valid)

    # print(confusion_matrix(y_valid,predictions))
    print(classification_report(y_valid, predictions))
    accuracy = (100 * accuracy_score(y_valid, predictions))

    return accuracy


def perform_RFE(x, y):
    x_train = np.delete(x, len(qm_descriptors), axis=1)
    scaler = StandardScaler(copy=True, with_mean=True, with_std=True)
    scaler.fit(x_train)
    x_train = scaler.transform(x_train)
    my_SVM_classifier = SVC(kernel="rbf")
    selector = RFE(estimator=my_SVM_classifier, n_features_to_select=1, step=1)
    selector.fit(x_train, y)

    return selector.ranking_


def build_SVM(filename, option, svm_type=None, poly_degree=None):
    # LOAD DATA
    descriptors = qm_descriptors
    X, Y = data_preprocess.load_data(filename, descriptors)

    if svm_type == None:
        svm_type = 'linear'

    if poly_degree == None:
        poly_degree = 2
        # print('training polynomial SVM of degree', poly_degree)

    if option == 'default':

        print('Start to train')
        print('*-----------------------------*')
        print('Start to train the default parameters.')

        accuracies_default = []
        for i in range(10):
            x_train, x_valid, y_train, y_valid = data_preprocess.split_data(X, Y, partition=0.20)
            accuracies_default.append(train_SVM(x_train, y_train, x_valid, y_valid))

        print('Average accuracy: %.2f' % np.mean(accuracies_default))


    elif option == 'train':

        print('*-----------------------------*')
        print('Find the best parameters.')

        params = []
        accuracies = []

        for i in range(10):
            x_train, x_valid, y_train, y_valid = data_preprocess.split_data(X, Y, partition=0.20)
            best_parameters = scan_parameters(x_train, y_train)
            params.append(best_parameters)
            accuracy = train_SVM(x_train, y_train, x_valid, y_valid, best_parameters)
            accuracies.append(accuracy)

        print('*-----------------------------*')
        print('Result.')
        print('*-----------------------------*')

        for i in range(len(accuracies)):
            print('Run ' + str(i + 1) + ' ', params[i], ' : ', accuracies[i])


    elif option == 'RFE':

        print('*-----------------------------*')
        print('Recursive feature estimation.')
        
        ranking = perform_RFE(X, Y)

        print('*-----------------------------*')
    
        print('*-----------------------------*')
        for d in range(len(qm_descriptors)):
            print(qm_descriptors[d], ranking[d])


    elif option == 'test':

        print('Testing......')
        print('*-----------------------------*')

        # kernels  = 'rbf'
        # Cs = 1
        # gammas = 1
        # degrees = 3
        # weights = None

        kernels = 'rbf'
        Cs = 10
        gammas = 0.1
        degrees = 3
        weights = None

        params_dict = {'kernel': kernels, 'C': Cs, 'class_weight': weights, 'degree': degrees, 'gamma': gammas}

        acc_list = []

        for i in range(10):
            x_train, x_valid, y_train, y_valid = data_preprocess.split_data(X, Y, partition=0.20)

            acc_list.append(train_SVM(x_train, y_train, x_valid, y_valid, params_dict))

        print('Result.')
        print('*-----------------------------*')
        print('Average accuracy: %.2f' % np.mean(acc_list))



    # execute only if run as a script
    import sys

    """
    svm types are :
        linear  [default]
        poly    [requires to specify the degree of polynomial]
        rbf     [gaussian]
        sigmoid [only for bi-classifier]

    """
option_choices = ["default","train","REF","test"]
for i in option_choices:
    build_SVM(r" ",i,{'kernel': ['linear', 'rbf', 'sigmoid']})
  
