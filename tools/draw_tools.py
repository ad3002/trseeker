#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 26.09.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
'''
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('Agg')

def draw_distribution_plot(distribution, image_file):
    ''' Draw distribution plot.
    '''
    x = []
    y = []
    for k, v in distribution.items():
        x.append(k)
        y.append(v) 
    x = np.asarray(x)
    y = np.asarray(y)
    plt.plot(x,y)
    print "Save image to", image_file
    plt.savefig(image_file)