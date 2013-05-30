#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2007-2009 Aleksey Komissarov ( ad3002@gmail.com )
# All rights reserved.
#
# This software is licensed as described in the file COPYING, which
# you should have received as part of this distribution.
"""
TRF search wrapper

- trf_search(file_name="")
- trf_search_in_dir(folder, verbose=False, file_suffix=".fa")

Command example: **wgs.AADD.1.gbff.fa 2 5 7 80 10 50 2000 -m -f -d -h**
"""

import os, shutil
from PyExp import sc_iter_filepath_folder
from trseeker.settings import load_settings

settings = load_settings()

def trf_search(file_name):
    """ TRF search in fasta file."""

    # TRF search, only data, mask
    params = "%s %s %s %s %s %s %s" % (settings["trf_settings"]["trf_match"],
                                       settings["trf_settings"]["trf_mismatch"],
                                       settings["trf_settings"]["trf_indel"],
                                       settings["trf_settings"]["trf_p"],
                                       settings["trf_settings"]["trf_q"],
                                       settings["trf_settings"]["trf_threshold"],
                                       settings["trf_settings"]["trf_length"],
                                       )
    if settings["trf_settings"]["trf_masked_file"]:
        params += " -m"
    if settings["trf_settings"]["trf_flanked_data"]:
        params += " -f"
    if settings["trf_settings"]["trf_data_file"]:
        params += " -d"
    if settings["trf_settings"]["trf_nohtml"]:
        params += " -h"
    if  settings["trseeker"]["os"] == "win":
        string = settings["trf_settings"]["trf_location"] + "trf400.dos.exe %s %s" % (file_name, params)
    else:
        if settings["trseeker"]["os"] == "WIN":
            string = settings["trf_settings"]["trf_location_WIN"] +" %s %s" % (file_name, params)
        else:
            string = settings["trf_settings"]["trf_location"] +" %s %s" % (file_name, params)
    print string
    os.system(string)

def trf_search_in_dir(folder, verbose=False, file_suffix=".fa", output_folder=None):
    """ TRF search in each file in the given folder.
    
    - folder with fasta files
    - verbose
    - file_suffix, function takes only files with suffix, default *.fa*
        
    """
    for file_name in sc_iter_filepath_folder(folder, mask=file_suffix):

        # check existence in output folder
        name = os.path.basename(file_name)
        if output_folder:
            if not os.path.isdir(output_folder):
                os.makedirs(output_folder)
            result_name = "%s.%s.%s" % (name, settings["trf_settings"]["trf_param_postfix"], "dat")
            dist_file = os.path.join(output_folder, result_name)


            result_file = os.path.abspath(result_name)
            if os.path.isfile(result_file):
                continue
            if os.path.isfile(dist_file):
                continue

        
        if verbose:
            print "Start trf search in %s" % file_name
        
        trf_search(file_name)
        if output_folder:
            dist_file = os.path.join(output_folder, result_name)
            print "Copy: ", result_file, dist_file
            shutil.copyfile(result_file, dist_file)
            os.unlink(result_file)

def _filter_by_bottom_array_length(obj, cutoff):
    if obj.trf_array_length > cutoff:
        return True
    else:
        return False

def _filter_by_bottom_unit_length(obj, cutoff):
    if obj.trf_period > cutoff:
        return True
    else:
        return False

def trf_filter_by_array_length(trf_file, output_file, cutoff):
    """ Create output TRF file with tandem repeats with length greater than from input file.
    Function returns number of tandem repeats in output file.
    """
    i = 0
    with open(output_file, "w") as fw:
        for obj in trf_reader(trf_file):
            if _filter_by_bottom_array_length(obj, cutoff):
                i += 1
                fw.write(obj.get_string_repr())
    print i
    return i

def trf_filter_by_monomer_length(trf_file, output_file, cutoff):
    """ Create output TRF file with tandem repeats with unit length greater than from input file.
    Function returns number of tandem repeats in output file.
    """
    i = 0
    with open(output_file, "w") as fw:
        for obj in trf_reader(trf_file):
            if _filter_by_bottom_unit_length(obj, cutoff):
                i += 1
                fw.write(obj.get_string_repr())
    print i
    return i

def trf_filter_exclude_by_gi_list(trf_file, output_file, gi_list_to_exclude):
    """ Create output TRF file with tandem repeats with GI that don't match GI_LIST
    List of GI, see TRF and FA specifications, GI is first value in TRF row.
    """
    with open(output_file, "w") as fw:
        for obj in trf_reader(trf_file):
            if not obj.trf_gi in gi_list_to_exclude:
                    fw.write(obj.get_string_repr())

def trf_representation(trf_file, trf_output, representation):
    """ Write TRF file tab delimited representation.
    representation: numerical|index|agc_apm|with_monomer|family
    """
    with open(trf_output, "w") as fw:
        for obj in trf_reader(trf_file):
            if representation == 'numerical':
                line = obj.get_numerical_repr()
            elif representation == 'index':
                line = obj.get_index_repr()
            elif representation == 'agc_apm':
                line = "%.2f\t%.2f\n" % (obj.trf_array_gc, obj.trf_pmatch)
            elif representation == "with_monomer":
                line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (obj.trf_id,
                                                     obj.trf_period,
                                                     obj.trf_array_length,
                                                     obj.trf_array_gc,
                                                     obj.trf_pmatch,
                                                     obj.trf_gi,
                                                     obj.trf_consensus)
            elif representation == 'family':
                line = obj.get_family_repr()

            fw.write(line)

def trf_write_field_n_data(trf_file, file_output, field, field_format="%s"):
    """ Write statistics data: field, N.
    """
    result = {}
    with open(file_output, "w") as fw:
        for obj in trf_reader(trf_file):
            value = field_format % getattr(obj, field)
            result.setdefault(value)
            result[value] += 1
        result = [ (value, n) for value, n in result.items() ]
        result.sort()
        for value, n in result:
            line = field_format + "\t%s\n"
            line = line % (value, n)
            fw.write(line)

def trf_write_two_field_data(trf_file, file_output, field_a, field_b):
    """ Write statistics data: field_a, field_b.
    """
    result = []
    with open(file_output, "w") as fw:
        for obj in trf_reader(trf_file):
            value_a = getattr(obj, field_a)
            value_b = getattr(obj, field_b)
            result.append([value_a, value_b])
        result.sort()
        for value_a, value_b in result:
            line = "%s\t%s\n" % (value_a, value_b)
            fw.write(line)

def count_trs_per_chrs(all_trf_file):
    """ Function prints chr, all trs, 3000 trs, 10000 trs
    """
    chr2n = {}
    chr2n_large = {}
    chr2n_xlarge = {}
    for trf_obj in trf_reader(all_trf_file):
        chr = trf_obj.trf_chr
        chr2n.setdefault(chr, 0)
        chr2n_large.setdefault(chr, 0)
        chr2n_xlarge.setdefault(chr, 0)
        chr2n[chr] += 1
        if trf_obj.trf_array_length > 3000:
            chr2n_large[chr] += 1
        if trf_obj.trf_array_length > 10000:
            chr2n_xlarge[chr] += 1
    for chr in chr2n:
        print chr, chr2n[chr], chr2n_large[chr], chr2n_xlarge[chr]

def count_trf_subset_by_head(trf_file, head_value):
    """ Function prints number of items with given fasta head fragment
    """
    total = 0
    n = 0
    total_length = 0
    for trf_obj in trf_reader(trf_file):
        total += 1
        if head_value in trf_obj.trf_head:
            n += 1
            total_length += trf_obj.trf_array_length
    return n, total, total_length

def fix_chr_names(trf_file, temp_file_name=None, case=None):
    """ Some fasta heads impossible to parse, so it is simpler to fix them postfactum
    
    Cases:
    
    - chromosome names in MSGC genome assembly [MSGC sequences]
    
    """

    if not temp_file_name:
        temp_file_name = trf_file + ".tmp"
    if case == "MSGC sequences":
        with open(temp_file_name, "a") as fw:
            for obj in trf_reader(trf_file):
                obj.trf_chr = obj.trf_gi
                fw.write(obj.get_string_repr())
    if os.path.isfile(temp_file_name):
        os.remove(trf_file)
        os.rename(temp_file_name, trf_file)