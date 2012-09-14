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

def cf_tab_trf(data):
    """ Transfrom tab data to trf_obj."""
    if len(data) > 21:
        # back compatibility case
        obj = TRModel()
        obj.trf_id = int(data[0])
        obj.trf_l_ind = int(data[1])
        obj.trf_r_ind = int(data[2])
        obj.trf_period = int(data[3])
        obj.trf_n_copy = float(data[4])
        obj.trf_l_cons = int(data[5])
        obj.trf_pmatch = int(data[6])
        obj.trf_indels = int(data[7])
        obj.trf_score = int(data[8])
        obj.trf_n_a = int(data[9])
        obj.trf_n_t = int(data[10])
        obj.trf_n_c = int(data[11])
        obj.trf_n_g = int(data[12])
        obj.trf_entropy = float(data[13])
        obj.trf_consensus = data[14]
        obj.trf_array = data[15]
        obj.trf_array_gc = float(data[16])
        obj.trf_consensus_gc = float(data[17])
        obj.trf_gi = data[18]
        obj.trf_head = data[19]
        obj.trf_param = data[20]
        obj.trf_array_length = int(data[21])
        obj.trf_chr = data[22]
        if len(data) > 23:
            obj.trf_joined = int(data[23])
        return obj
    else:
        obj = TRModel()
        obj.trf_id = int(data[0])
        obj.trf_l_ind = int(data[1])
        obj.trf_r_ind = int(data[2])
        obj.trf_period = int(data[3])
        obj.trf_n_copy = float(data[4])
        obj.trf_pmatch = int(data[5])
        obj.trf_consensus = data[6]
        obj.trf_array = data[7]
        obj.trf_array_gc = float(data[8])
        obj.trf_consensus_gc = float(data[9])
        obj.trf_gi = data[10]
        obj.trf_head = data[11]
        obj.trf_param = data[12]
        obj.trf_array_length = int(data[13])
        obj.trf_chr = data[14]
        obj.trf_joined = int(data[15])
        return obj

def trf_to_array_fasta(trf_file, fa_output):
    """ Write TRF file like tandem repeat array fasta."""
    reader = TabDelimitedFileIO(format_func=cf_tab_trf)
    with open(fa_output, "w") as fw:
        for trf_obj in reader.read_online(trf_file):
            fw.write(trf_obj.get_fasta_repr())

def trf_to_monomer_fasta(trf_file, fa_output):
    """ Write TRF file like tandem repeat array fasta."""
    reader = TabDelimitedFileIO(format_func=cf_tab_trf)
    with open(fa_output, "w") as fw:
        for trf_obj in reader.read_online(trf_file):
            fw.write(trf_obj.get_monomer_fasta_repr())

def get_trid2trf_obj(trf_large_file):
    ''' Read TRs file and return trid->trf_obj dictionary.'''
    trid2obj = {}
    for trf_obj in sc_iter_tab_file(trf_large_file, TRModel):
        trid2obj[trf_obj.trf_id] = trf_obj
    return trid2obj


def read_trid2meta(file_name):
    ''' Load trid to full index dictionary.'''

    print "Load trid to full index dictionary"
    trid2meta = {}
    for trf_obj in sc_iter_tab_file(file_name, TRModel):
        trid2meta[trf_obj.trf_id] = str(trf_obj)
    return trid2meta

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

def get_all_trf_objs(trf_large_file):
    """ Return list of trf_obj from given trf_large_file."""

    result = []
    for trf_obj in sc_iter_tab_file(trf_large_file, TRModel):
        result.append(trf_obj)
    return result

def get_trfid_obj_dict(trf_large_file):
     """ Return dicionary trf_id totrf_obj from given trf_large_file."""
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

