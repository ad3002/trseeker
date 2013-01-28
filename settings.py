#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 06.09.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
Trseeker settings loader.
'''

import yaml, os

SETTINGS_FILENAME = "settings.yaml"
NGRAM_LENGTH = 21
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
