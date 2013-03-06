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

class KmerSliceModel(AbstractModel):
    '''
    '''
    dumpable_attributes = ["kmer",
                             "local_tf",
                             "df",
                             ]
    int_attributes = ["local_tf",
                             "df",
                        ]

    def __init__(self, data):
      super(KmerSliceModel, self).__init__()
      a, b, c = data.split(":")
      self.kmer = a
      self.local_tf = int(b)
      self.df = int(c)

    def set_with_dict(self, dictionary):
      raise NotImplemented 

    def set_with_list(self, dictionary):
      raise NotImplemented

    def __str__(self):
      return "%s:%s:%s" % (self.kmer, self.local_tf, self.df)


class SliceTreeModel(AbstractModel):
    ''' Model for kmers tree node.
    '''

    dumpable_attributes = [
                           "deep",
                           "size",
                           "blast_fams",
                           "maxdf",
                           "nmaxdf",
                           "pmaxdf",
                           "gc_var",
                           "units",
                           "trs",
                           "kmers",
                           ]

    int_attributes = [
                           "deep",
                           "size",
                           "maxdf",
                           "nmaxdf",
                      ]

    float_attributes = [
                        "pmaxdf",
                        "gc_var",
    ]

    list_attributes = ["blast_fams",
                       "units",
                       "trs",
                       "kmers",
                       ]

    list_attributes_types = {"blast_fams":str,
                             "units":int,
                             "trs":int,
                             "kmers":KmerSliceModel,
                             }

    

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