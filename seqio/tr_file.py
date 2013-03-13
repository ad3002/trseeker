#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 07.09.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
Working with tab-delimited TRs datasets files.
TODO: check it.
'''
import os
from trseeker.seqio.tab_file import TabDelimitedFileIO, sc_iter_tab_file, \
    sc_iter_simple_tab_file
from trseeker.models.trf_model import TRModel

def read_trid2ngrams(annotation_ngram_folder, trf_large_file):
    """ Read trid -> [(ngram, tf), (rev_ngram, tf), ...] data."""

    trid2ngrams = {}
    for trf_obj in sc_iter_tab_file(trf_large_file, TRModel):

        file_name = os.path.join(annotation_ngram_folder, "%s.ngram" % trf_obj.trf_id)
        trid2ngrams.setdefault(trf_obj.trf_id, [])

        for data in sc_iter_simple_tab_file(file_name):
            ngram = data[0]
            rev_ngram = data[1]
            tf = float(data[2])
            trid2ngrams[trf_obj.trf_id].append((ngram, tf))
            trid2ngrams[trf_obj.trf_id].append((rev_ngram, tf))
    return trid2ngrams

def read_trid2meta(file_name):
    ''' Load trid to full index dictionary as string.'''

    print "Load trid to full index dictionary"
    trid2meta = {}
    for trf_obj in sc_iter_tab_file(file_name, TRModel):
        trid2meta[trf_obj.trf_id] = str(trf_obj)
    return trid2meta

def get_all_trf_objs(trf_large_file):
    """ Return list of trf_obj from given trf_large_file."""
    result = []
    for trf_obj in sc_iter_tab_file(trf_large_file, TRModel):
        result.append(trf_obj)
    return result

def get_trf_objs_dict(trf_large_file):
    """ Return dict of trf_obj from given trf_large_file."""
    result = {}
    for trf_obj in sc_iter_tab_file(trf_large_file, TRModel):
        result[trf_obj.trf_id] = trf_obj
    return result

def get_trfid_obj_dict(trf_large_file):
     """ Return dicionary trf_id to trf_obj from given trf_large_file."""
     trs_dataset = get_all_trf_objs(trf_large_file)
     trid2obj = {}
     for trf_obj in trs_dataset:
        trid2obj[trf_obj.trf_id] = trf_obj
     return trid2obj

def save_trs_dataset(trs_dataset, output_file):
    """ Save trs dataset to file."""
    with open(output_file, "w") as fh:
        for trf_obj in trs_dataset:
            data = str(trf_obj)
            fh.write(data)

def save_trs_as_fasta(trf_file, fasta_file):
    ''' Save TRs dataset as one fasta file.
    '''
    trf_objs = []
    for trf_obj in sc_iter_tab_file(trf_file, TRModel):
        trf_objs.append(trf_obj)
    with open(fasta_file, "w") as fh_fasta:
        for trf_obj in trf_objs:
            fh_fasta.write(trf_obj.get_fasta_repr())

