#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 26.09.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
Simple statistics.
'''
import math
from collections import defaultdict

def get_variance(data):
    '''Calculated variance for given list.

    ( sum(X*X) - sum(X) * mean_x ) / (n-1)
    '''
    n = 0
    Sum = 0
    Sum_sqr = 0

    if not data:
        raise "Empty arrays for variance computation."
    n = len(data)
    if n == 1:
        return 0
    sum_x = sum(data)
    sum_xx = sum([x*x for x in data])
    mean = float(sum_x) / n
    variance = (sum_xx - sum_x * mean) / (n - 1)
    return variance

def get_sigma(data):
    ''' Calculate sigma, return sum(module(xi - mean).
    '''
    if not data:
        raise "Empty arrays for variance computation."
    n = len(data)
    if n == 1:
        return 0
    mean = get_mean(data)
    return sum( [ abs(x - mean) for x in data ] )

def get_mean(data):
    ''' Calculated mean for given list.
    '''
    if not data:
        raise Exception("Empty data.")
    sum_x = sum(data)
    mean = float(sum_x) / len(data)
    return mean

def get_standard_deviation(variance):
    ''' Get sample deviation for variance.
    '''
    if variance<0:
        raise "Wrong variance value %s" % variance
    return math.sqrt(variance)

def t_test(sample_mean, dist_mean, variance, N):
    ''' T test.
    
    - sample_mean: sample mean
    - dist_mean: distribution mean
    - variance: sample variance
    - N: sample size
    '''
    if N<=0:
        raise "Wrong N value %s" % N
    if variance<=0:
        raise "Wrong variance value %s" % variance
    return (sample_mean-dist_mean)/float( math.sqrt( variance/N ) )

def get_element_frequences(data):
    ''' Get defaultdictionary of elements frequences in given list
    TODO: refr this.
    '''
    d = defaultdict(int)
    for element in data:
        d[element] += 1
    return d

def get_simple_statistics(data):
    ''' Return dictionary containing simple statistics for given list. 
    Dictionary keys: mean, variance, sigma, sample_derivation.
    
    >>> result = {
    >>>    'mean': get_mean(data),
    >>>    'variance': get_variance(data),
    >>>    'sigma': get_sigma(data),
    >>>    'standard_deviation': get_standard_deviation(data),
    >>> }
    
    '''
    variance = get_variance(data)
    mean = get_mean(data)
    result = {
        'mean': mean,
        'variance': variance,
        'sigma': get_sigma(data),
        'standard_deviation': get_standard_deviation(variance),
    }
    return result