#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 23.02.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
Functions related to ngram (k-mer).
'''
from collections import defaultdict

def generate_ngrams(text, n=12):
    ''' Yields all ngrams of length n from given text. 
    
    - n: ngram length
    '''
    for i in xrange(0, len(text) - n + 1):
        yield i, text[i:i + n]

def generate_window(text, n=12, step=None):
    ''' Yields all ngrams of length n from given text. 
    - n: ngram length
    - step: window step
    '''
    if not step:
        step = n / 2
    for i in xrange(0, len(text) - n + 1, step):
        yield i, text[i:i + n]


def get_ngrams(text, m=5, n=12):
    ''' Returns m most frequent (ngram of length n, tf) tuples for given text.
    
    - m: number of returned ngrams
    - n: ngram length
    '''

    ngrams = {}
    for pos, ngram in generate_ngrams(text, n=n):
        ngrams.setdefault(ngram, 0)
        ngrams[ngram] += 1
    ngrams = [(key, value) for key, value in ngrams.items()]
    ngrams.sort(reverse=True, key=lambda x: x[1])
    return ngrams[:m]

def get_ngrams_freq(text, m=5, n=12):
    ''' Returns m most frequent (ngram of length n, fraction of possible ngrams) tuples for given text.
    
    - m: number of returned ngrams
    - n: ngram length
    '''

    ngrams = defaultdict(int)
    for pos, ngram in generate_ngrams(text, n=n):
        ngrams[ngram] += 1
    text_length = float(len(text) - n + 1)
    ngrams = [(key, value, value / text_length) for key, value in ngrams.items()]
    ngrams.sort(reverse=True, key=lambda x: x[1])
    return ngrams[:m]

def get_ngrams_feature_set(text, m=5, n=12):
    ''' Returns a feature set {'ngram':'ngram',...}  of m most frequent ngram of length n for given text.
    
    - m: number of returned ngrams
    - n: ngram length
    '''

    ngrams = {}
    for pos, ngram in generate_ngrams(text, n=n):
        ngrams.setdefault(ngram, 0)
        ngrams[ngram] += 1
    result = [(key, value) for key, value in ngrams.items()]
    ngrams.sort(reverse=True, key=lambda x: x[1])
    data = {}
    for key, value in enumerate(result[:m]):
        data[key] = key
    return data

def get_ngram_freq_distance(ngrams_a, ngrams_b):
    ''' Returns a distance between two ngram sets where distance is a sum(min(ngram_a, ngram_b) for each common ngram)
    
    - ngrams_a: dictionary {ngram:n, ...}
    - ngrams_b: dictionary {ngram:n, ...}
    '''
    distance = 0
    common_ngrams = [ngram for ngram, fr in ngrams_a.items() if ngram in ngrams_b]
    for ngram in common_ngrams:
        distance += min(ngrams_a[ngram], ngrams_b[ngram])
    return distance

def get_ngram_common_distance(ngrams_a, ngrams_b):
    ''' Returns a distance between two ngram sets where distance is a len()
    
    - ngrams_a: dictionary {ngram:n, ...}
    - ngrams_b: dictionary {ngram:n, ...}
    '''
    common_ngrams = [ngram for ngram, fr in ngrams_a.items() if ngram in ngrams_b]
    return len(common_ngrams)

