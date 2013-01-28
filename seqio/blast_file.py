#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 05.06.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com
'''
Working with BLAST tab-delimited output files.
'''
import os
from trseeker.models.blast_model import read_blast_file

def _get_blast_result_intervals(blast_file, length):
    ''' Return gi->sorted list of blast objects. '''

    # read data and Filter A: min_align
    gi_to_results = read_blast_file(blast_file, length)

    # create correct start/end positions
    for gi in gi_to_results:
        for i, blast_obj in enumerate(gi_to_results[gi]):
            if blast_obj.query_start > blast_obj.query_end:
                temp = blast_obj.query_start
                gi_to_results[gi][i].query_start = blast_obj.query_end
                gi_to_results[gi][i].query_end = temp

    # sort data
    for gi in gi_to_results:
        gi_to_results[gi].sort(key=lambda x: (x.query_start, x.query_end))

    return gi_to_results

def _remove_nested_join_overlap(blast_dataset):
    ''' Remove nested matches. 
        A: skip inners
            # -------
            #  ----
            # 
            # -------
            # ----
            # 
            # -------
            #    ---- 
            # 
            # -------
            # -------
        B: join overlapped
            # -------      -----
            #   ---------       -----
            
    '''
    for gi in blast_dataset:
        previous_item = None
        for i, item in enumerate(blast_dataset[gi]):
            if item is None:
                continue
            # check data
            assert item.query_start < item.query_end
            if not previous_item:
                previous_item = item
                pi = i
                continue
            # A: skip inners
            if item.query_start >= previous_item.query_start \
                    and item.query_end <= previous_item.query_end:
                blast_dataset[gi][i] = None
                continue

            # B: join overlapped
            if item.query_start <= previous_item.query_end \
                    and item.query_end > previous_item.query_end:

                blast_dataset[gi][pi].query_end = item.query_end
                previous_item.query_end = item.query_end

                blast_dataset[gi][pi].changed = True
                blast_dataset[gi][i] = None
                continue
            previous_item = item
            pi = i
        blast_dataset[gi] = [x for x in blast_dataset[gi] if x]
    return blast_dataset

def _join_gapped(blast_dataset, gap_size, min_align,):
    ''' Join gapped matches.
        C: remove short alignments 
        D: join gapped
        #  ------ 
        #            ------
    '''

    for gi in blast_dataset:
        previous_item = None
        for i, item in enumerate(blast_dataset[gi]):
            if item is None:
                continue

            # C: remove short alignments 
            length = item.query_end - item.query_start
            if length < min_align:
                blast_dataset[gi][i] = None
                continue

            if not previous_item:
                previous_item = item
                pi = i
                continue

            # D: join gapped
            if abs(previous_item.query_end - item.query_start) <= gap_size:
                blast_dataset[gi][pi].query_end = item.query_end
                previous_item.query_end = item.query_end
                blast_dataset[gi][pi].changed = True
                blast_dataset[gi][pi].gapped = True
                blast_dataset[gi][i] = None
                continue

            previous_item = item
            pi = i

        blast_dataset[gi] = [x for x in blast_dataset[gi] if x]
    return blast_dataset

def _filter_blast_dataset(blast_dataset, min_length):
    ''' Filter by minimal alignment length.

    - E: skip short aligns
    '''
    # E: skip short aligns 
    for gi in blast_dataset:
        for i, item in enumerate(blast_dataset[gi]):
            length = item.query_end - item.query_start
            if length < min_length:
                blast_dataset[gi][i] = None
                continue
        blast_dataset[gi] = [x for x in blast_dataset[gi] if x]
    blast_dataset = dict([(k, v) for k, v in blast_dataset.items() if v])
    return blast_dataset

def _format_output(blast_dataset, format_function):
    ''' Format blast_dataset output.'''
    assert hasattr(format_function, "__call__")
    return format_function(blast_dataset)

def format_function_repbase(blast_dataset):
    ''' Join with keys with ;.'''
    return ";".join(blast_dataset.keys())

def format_function_self(blast_dataset):
    ''' Join keys with ,.'''
    return ",".join(blast_dataset.keys())

def get_blast_result(blast_file, length, gap_size=1000, min_align=500, min_length=2400, format_function=None):
    """ Get blast result."""

    blast_dataset = _get_blast_result_intervals(blast_file, length)
    blast_dataset = _remove_nested_join_overlap(blast_dataset)
    blast_dataset = _join_gapped(blast_dataset, gap_size, min_align)
    blast_dataset = _filter_blast_dataset(blast_dataset, min_length)
    if not format_function:
        format_function = format_function_self
    result = _format_output(blast_dataset, format_function)
    return result

def update_with_repbase_blast_result(trs_dataset, annotation_self_folder, filters):
    ''' Add Repbase blast result to trs_dataset.'''
    for i, trf_obj in enumerate(trs_dataset):

        blast_output_file = os.path.join(annotation_self_folder, "%s.blast" % trf_obj.trf_id)
        result = get_blast_result(blast_output_file,
                                  trf_obj.trf_array_length,
                                  gap_size=filters["blast_gap_size"],
                                  min_align=filters["min_align"],
                                  min_length=filters["min_length"],
                                  format_function=format_function_repbase)
        trs_dataset[i].trf_repbase = result

    return trs_dataset

def update_with_self_blast_result(trs_dataset, annotation_self_folder, filters):
    ''' Add vs_self blast result to trs_dataset.'''
    n = len(trs_dataset)
    for i, trf_obj in enumerate(trs_dataset):
        print i, n, "\r",
        blast_output_file = os.path.join(annotation_self_folder, "%s.blast" % trf_obj.trf_id)
        with open(blast_output_file) as fh:
            data = fh.read()
        if data.startswith("ALPHA"):
            data = data.split()
            ref = int(data[1])
            trs_dataset[i].trf_family_self = ('ALPHA', ref)
            continue
        result = get_blast_result(blast_output_file,
                                  trf_obj.trf_array_length,
                                  gap_size=filters["blast_gap_size"],
                                  min_align=filters["min_align"],
                                  min_length=filters["min_length"],
                                  format_function=format_function_self)
        trs_dataset[i].trf_family_self = result
    print
    return trs_dataset

def update_with_ref_blast_result(trs_dataset, annotation_self_folder, filters):
    ''' Add vs_reference blast result to trs_dataset.'''
    n = len(trs_dataset)
    for i, trf_obj in enumerate(trs_dataset):
        print i, n, "\r",

        blast_output_file = os.path.join(annotation_self_folder, "%s.blast" % trf_obj.trf_id)
        result = get_blast_result(blast_output_file,
                                  trf_obj.trf_array_length,
                                  gap_size=filters["blast_gap_size"],
                                  min_align=filters["min_align"],
                                  min_length=filters["min_length"],
                                  format_function=format_function_self)
        trs_dataset[i].trf_family_ref = result
    print
    return trs_dataset

