# !/usr/bin/env python
# -*-coding:utf-8 -*-

# File       : file_path.py
# Description：This file is used for RF model training

import numpy
import pickle
import Template


def RF_model(filename, X_data):
    # LOAD DATA
    x = numpy.delete(X_data, len(Template.ALL_PROPERTIES), axis=1)

    print('Start')
    print('*-----------------------------*')

    my_forest = pickle.load(open(Template.RANDOM_FOREST_MODEL, 'rb'))
    print('Loaded Random Forest Classifier from:', Templatet.RANDOM_FOREST_MODEL)

    m = 'Results. \n'
    m = m + '*-----------------------------*\n'

    predicted_y = (my_forest.predict(x))
    class_probability = (my_forest.predict_proba(x))

    # print(str(X_data[0]))
    m = m + 'C_IDX,  PRED_ACTIVITY, PROBAB_of_INACT, PROBAB_of_ACT \n'
    for i in range(len(predicted_y)):
        m_i = ', '.join([str(int(X_data[i][-1])), str(predicted_y[i]), "{0:.2f}".format(class_probability[i, 0]),
                         "{0:.2f}".format(class_probability[i, 1])])
        m = m + m_i + '\n'

    # print(m)

    fo = open(filename, 'w')
    fo.write(m)
    fo.close()
    print('>> Results written to result file.')
    return predicted_y, class_probability



