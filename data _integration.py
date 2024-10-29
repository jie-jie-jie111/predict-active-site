# !/usr/bin/env python
# -*-coding:utf-8 -*-

# File       : data_integration.py
# Description：This file is used to consolidate all descriptors
# Output data format： name, label , Q, CHB, BO, ASA, Fukui, AN
# Output result: y：0——inactive，1——active


import csv
import os

import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split

qm_descriptors = ['Q', 'CHB', 'BO', 'ASA', 'Fukui', 'AN']


def load_data(data_file, descriptors):
    my_descriptors = qm_descriptors
    # print(len(my_descriptors))
    descriptors_data = {'labels': [], 'idx': []}
    for i in my_descriptors:
        descriptors_data[i] = []
    with open(data_file) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        l = 0
        for row in readCSV:
            descriptors_data['idx'].append(l)
            descriptors_data['Q'].append(float(row[2]))
            descriptors_data['CHB'].append(float(row[3]))
            descriptors_data['BO'].append(float(row[4]))
            descriptors_data['ASA'].append(float(row[5]))
            descriptors_data['Fukui'].append(float(row[6]))
            descriptors_data['AN'].append(float(row[7]))
            # descriptors_data['labels'].append(float(row[8]))
            l = l + 1

    # add descriptors columns to X in same order as in all_descriptors
    labels = np.array(descriptors_data['labels'])
    print(labels)
    X = np.zeros((len(my_descriptors) + 1, len(labels)))
    for i in my_descriptors:
        # print (i, my_descriptors.index(i))
        X[my_descriptors.index(i):] = descriptors_data[i]
    X[(len(my_descriptors)):] = descriptors_data['idx']
    X = X.transpose()
    # print(X[0])
    print('Reactivity site: %.2f' % (100 * sum(labels) / len(labels)))
    print('The dimension of X: ', (X.shape))
    print('The dimension of Y: ', (labels.shape))
    return X, labels


def do_down_sampling(x, y):
    # Mark 0 and 1 for each possible locus
    negatives = np.where(y == 0)[0]
    positives = np.where(y == 1)[0]

    n_negatives = len(negatives)
    n_positives = len(positives)
    print('!!Find', n_positives, 'positives and ', n_negatives, 'negatives.')
    print('Increase randomness and reduce duplication')
    negatives_downsampled = np.random.choice(negatives, size=n_positives, replace=False)
    print(len(negatives_downsampled))
    y_resampled = np.concatenate((y[positives], y[negatives_downsampled]))
    x_resampled = np.concatenate((x[positives], x[negatives_downsampled]))
    print('The dimension of new samples:', len(y_resampled))
    return x_resampled, y_resampled


def split_negatives_positives(x, label):
    negatives = np.where(label == 0)[0]
    positives = np.where(label == 1)[0]
    positives_x = np.zeros(len(positives))
    j = 0
    for index in positives:
        positives_x[j] = x[index]
        j = j + 1
    negatives_x = np.zeros(len(negatives))
    j = 0
    for index in negatives:
        negatives_x[j] = x[index]
        j = j + 1
    return positives_x, negatives_x


def split_data(X, Y, partition):
    x1, x2, y1, y2 = train_test_split(X, Y, test_size=partition)
    print('SET 1: \n\t dataset num: ' + str(len(y1)) + '\n\t probability: %.2f' % (100 * sum(y1) / len(y1)))
    print('SET 2: \n\t dataset num: ' + str(len(y2)) + '\n\t probability: %.2f' % (100 * sum(y2) / len(y2)))
    return x1, x2, y1, y2



