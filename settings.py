#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 06.09.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
Trseeker settings loader.

Raw parameters:

>>> SETTINGS_FILENAME = "settings.yaml"
>>> NGRAM_LENGTH = 12
>>> NGRAM_N = 100000

Yaml file:

>>> trseeker:
>>>     os: linux64
>>>     root_dir: /root/Dropbox/workspace/trseeker
>>>     work_dir: /home
>>> blast_settings:
>>>     blast_location: 
>>>     repbase_db_folder: /home/rebase_blast_db/repbase
>>>     blast_e: 0.000000000001
>>>     blast_b: 1000000
>>>     blast_v: 1000000
>>> trf_settings:
>>>     trf_location: /root/trf404.linux64
>>>     trf_match: 2
>>>     trf_mismatch: 5
>>>     trf_indel: 7
>>>     trf_p: 80
>>>     trf_q: 10
>>>     trf_threshold: 50
>>>     trf_length: 2000
>>>     trf_param_postfix: 2.5.7.80.10.50.2000
>>>     trf_masked_file: False
>>>     trf_flanked_data: True
>>>     trf_data_file: True
>>>     trf_nohtml: True
>>>     overlapping_cutoff: 10
>>>     overlapping_cutoff_proc: 30
>>>     overlapping_gc_diff: 0.05
>>> ngrams_settings:
>>>     ngram_length: 12
>>>     ngram_m: 100000 

'''

import yaml, os

SETTINGS_FILENAME = "settings.yaml"

NGRAM_LENGTH = 12
NGRAM_N = 100000000

def load_settings():
    ''' Load settings from yaml file.'''
    file = os.path.abspath(__file__)
    settings_file = os.path.join(os.path.split(file)[0],
                                 SETTINGS_FILENAME)
    with open(settings_file) as fh:
        settings = yaml.load(fh)
    return settings

def save_settings(settings):
    ''' Save settings to yaml file.'''
    file = os.path.abspath(__file__)
    settings_file = os.path.join(os.path.split(file)[0],
                                 SETTINGS_FILENAME)
    with open(settings_file, "w") as fh:
        yaml.dump(fh, settings)
