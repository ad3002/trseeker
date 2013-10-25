#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 07.06.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
"""
Function related to TRs classification.
"""

from trseeker.seqio.tab_file import sc_iter_tab_file
from trseeker.tools.sequence_tools import fix_strand
from trseeker.models.trf_model import TRModel
from trseeker.models.trf_model import TRsClassificationModel


def get_trs_types(trf_all_file, trf_all_class_file, settings):
    """
    Classify TRs into raw types.
    @param trf_all_file: input file with TRs
    @param trf_all_class_file: input file with TRs classification
    @param settings: PySatDNA settings dictionary
    @return: (trs_types, class_objs, trf_objs)
    """
    trs_types = {
        "all": 0,
        "micro": 0,
        "perfect": 0,
        "too_short": 0,
        "100bp": 0,
        "gc": 0,
        "x4": 0,
        "entropy": 0,
        "ok": 0,
    }

    print "Load classification"
    class_objs = []
    for class_obj in sc_iter_tab_file(trf_all_class_file, TRsClassificationModel):
        class_objs.append(class_obj)
    
    print "Get   types"
    trf_objs = []
    for i, trf_obj in enumerate(sc_iter_tab_file(trf_all_file, TRModel)):
        fixed_array = fix_strand(trf_obj.trf_array)
        if trf_obj.trf_array != fixed_array:
            trf_obj.trf_array = fixed_array
            trf_obj.trf_consensus = fix_strand(trf_obj.trf_consensus)
        trs_types["all"] += 1
        ok = True
        if trf_obj.trf_array_length < settings["other"]["trs_types"]["possible_array_length"]:
            trs_types["too_short"] += 1
            class_objs[i].class_100bp = "100bp"
            ok = False
        if len(trf_obj.trf_consensus) < settings["other"]["trs_types"]["microsat_monomer"]:
            trs_types["micro"] += 1
            class_objs[i].class_micro = "micro"
            ok = False
        if trf_obj.trf_pmatch == 100.0:
            trs_types["perfect"] += 1
            class_objs[i].class_perfect = "perfect"
            ok = False
        if trf_obj.trf_array_length < settings["other"]["trs_types"]["min_array_length"]:
            trs_types["100bp"] += 1
            class_objs[i].class_100bp = "100bp"
            ok = False
        if trf_obj.trf_array_gc < settings["other"]["trs_types"]["min_gc"] or trf_obj.trf_array_gc > settings["other"]["trs_types"]["max_gc"]:
            trs_types["gc"] += 1
            class_objs[i].class_gc = "gc"
            ok = False
        if trf_obj.trf_n_copy < settings["other"]["trs_types"]["min_copy_number"]:
            trs_types["x4"] += 1
            class_objs[i].class_x4 = "x4"
            ok = False
        if trf_obj.trf_entropy < settings["other"]["trs_types"]["min_entropy"]:
            trs_types["entropy"] += 1
            class_objs[i].class_entropy = "entropy"
            ok = False
        if ok:
            trs_types["ok"] += 1
            class_objs[i].class_good = "good"
        print trs_types["all"], trs_types["micro"], trs_types["perfect"], trs_types["100bp"], trs_types["gc"], trs_types["x4"], trs_types["entropy"], trs_types["ok"], "\r",         
        trf_objs.append(trf_obj)
    print
    return trs_types, class_objs, trf_objs