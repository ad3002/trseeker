#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 30.06.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

from PyExp.models.abstract_model import AbstractModel
from trseeker.seqio.tab_file import sc_iter_tab_file, sc_iter_simple_tab_file

class NgramModel(AbstractModel):
    ''' Ngram model.

    Dumpable attributes:

    - id (int)
    - rev_id (int)
    - ngram
    - rev_ngram
    - tf (int), ngram frequency
    - df (int), document frequency
    - edf (int), ngram frequency for both strands
    - etf (int)
    '''

    dumpable_attributes = ["id",
                           "rev_id",
                           "ngram",
                           "rev_ngram",
                           "tf",
                           "df",
                           "etf",
                           "edf",
                           ]
    int_attributes = ["id",
                      "rev_id",                      
                      "df",
                      "edf",
                      ]

    float_attributes = ["tf","etf"]

class NgramToTRModel(AbstractModel):
    ''' Model for ngram to TR data.

    Dumpable attributes:

    - id (int)
    - trids (list of int)
    - tfs (list of int)

    Other attributes:

    - idtf (list)

    '''
    dumpable_attributes = ["id",
                           "trids",
                           "tfs",
                           ]
    int_attributes = ["id",
                      ]

    float_attributes = []

    list_attributes = ["trids",
                       "tfs",
                       ]

    list_attributes_types = {"trids":int,
                             "tfs":int,
                             }

    other_attributes = {"idtf":[]}

def sc_ngram_reader(file_name):
    ''' Read file with ngrams data.'''
    for obj in sc_iter_tab_file(file_name, NgramModel):
        yield obj

def sc_ngram_trid_reader(file_name):
    ''' Read file with ngram to TRs+tfs data.'''
    for nid, sid_list, tf_list  in sc_iter_simple_tab_file(file_name):
        sids = sid_list.split(",")
        tfs = tf_list.split(",")

        result = []
        for i in xrange(0, len(sids)):
            result.append((int(sids[i]), float(tfs[i])))
        yield int(nid), result
