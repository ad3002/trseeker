#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 07.09.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com
"""
Jellyfish python wrapper

Used settings:
settings["blast_settings"]["jellyfish_location"]
settings["jellyfish_settings"]["hash_size"]
settings["jellyfish_settings"]["threads"]
settings["jellyfish_settings"]["both_strands"]
"""

import os
from trseeker.settings import load_settings
from PyExp import sc_iter_filepath_folder
import subprocess
from trseeker.tools.seqfile import sort_file_by_int_field
from trseeker.tools.ngrams_tools import process_list_to_kmer_index
from trseeker.tools.sequence_tools import get_revcomp

settings = load_settings()
location = settings["blast_settings"]["jellyfish_location"]


def count_kmers(input_file, output_prefix, k, mintf=None):
    """
    Count kmers with Jellyfish.
    @param input_file: input fasta or fastq file name
    @param output_prefix: output path with file prefix
    @param k: k-mer length
    @param mintf: count only kmers with frequency greater than mintf
    @return: None
    """
    params = {
        "location": location,
        "input_fasta": input_file,
        "k": k,
        "hash_size": settings["jellyfish_settings"]["hash_size"],
        "hash_bits": settings["jellyfish_settings"]["hash_bits"],
        "threads": settings["jellyfish_settings"]["threads"],
        "both_strands": settings["jellyfish_settings"]["both_strands"],
        "output_prefix": output_prefix,
        "mintf": "",
    }
    if mintf:
        params["mintf"] = "--lower-count=%s" % mintf
    command = "%(location)s count %(mintf)s -m %(k)s -o %(output_prefix)s -c %(hash_bits)s -s %(hash_size)s %(both_strands)s -t %(threads)s %(input_fasta)s" % params
    print "Execute:", command
    os.system(command)


def merge_kmers(folder, output_prefix, output_file):
    """
    Merge Jellyfish count output to one file.
    @param folder: folder with input files
    @param output_prefix: output prefix
    @param output_file: output file
    @return: None
    """
    if not output_file.endswith(".jf"):
        output_file += ".jf"
    params = {
        "location": location,
        "output_file": output_file,
        "output_prefix": output_prefix,
    }
    file_count = 0
    for file_name in sc_iter_filepath_folder(folder, mask="."):
        if output_prefix in file_name:
            file_count += 1
    assert file_count > 0
    if file_count == 1:
        command = "cp %(output_prefix)s_0 %(output_file)s" % params
    else:
        command = "%(location)s merge -o %(output_file)s %(output_prefix)s\_*" % params
    print "Execute:", command
    os.system(command)
    command = "Execute: rm %(output_prefix)s_*" % params
    print command
    os.system(command)


def stats_kmers(db_file, stats_file):
    """
    Compute statistics for kmers.
    @param db_file: jf db file
    @param stats_file: output stats file
    @return: None
    """
    params = {
        "location": location,
        "db_file": db_file,
        "stats_file": stats_file,
    }
    command = "%(location)s stats --verbose -o %(stats_file)s %(db_file)s" % params
    print "Execute:", command
    os.system(command)


def histo_kmers(db_file, histo_file):
    """
    Compute frequencies histogram.
    @param db_file: jf db file
    @param histo_file: histogram output file
    @return: None
    """
    params = {
        "location": location,
        "db_file": db_file,
        "histo_file": histo_file,
    }
    command = "%(location)s histo --verbose -o %(histo_file)s %(db_file)s" % params
    print "Execute:", command
    os.system(command)


def dump_kmers(db_file, kmers_file, dumpmintf):
    """
    Dump and sort k-mers database to tab-delimited file.
    @param db_file: jf db file
    @param kmers_file: output kmer file
    @param dumpmintf: minimum tf to dump
    @return: None
    """
    params = {
        "location": location,
        "db_file": db_file,
        "kmers_file": kmers_file,
        "dumpmintf": dumpmintf,
    }
    command = "%(location)s dump --column -L %(dumpmintf)s --tab -o %(kmers_file)s %(db_file)s" % params
    print "Execute:", command
    os.system(command)
    sort_file_by_int_field(kmers_file, 1)


def query_kmers(db_file, query_hashes, both_strands=True, verbose=True):
    """
    Query jellyfish database.
    @param db_file: jf db file
    @param query_hashes: kmers to query
    @param both_strands: use both strands
    @param verbose: verbose
    @return: dictionary hash to tf
    """
    params = {
        "location": location,
        "db_file": db_file,
        "query_hashes": query_hashes,
        "both_strands": "",
    }
    if both_strands:
        params["both_strands"] = "-C"
    command = "%(location)s query %(both_strands)s %(db_file)s" % params
    if verbose:
        print command
    final_result = {}
    n = len(query_hashes)
    for k in xrange(0,n,100):
        pp = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE, universal_newlines=True)
        for i, query in enumerate(query_hashes[k:k+100]):
            pp.stdin.write(query+" ")
        data = pp.communicate()
        error = data[1]
        if "Can't open file" in error:
            return None
        for item in data[0].strip().split("\n"):
            if item:
                if verbose:
                    print item
                key, value = item.strip().split()
                final_result[key] = value
            else:
                if verbose:
                    print data
                final_result[-1] = data[1]
    return final_result


def query_and_write_coverage_histogram(db_file, query_sequence, output_file, k=23):
    """
    Save coverage histogram into output_file for given query_sequence.
    @param db_file:
    @param query_sequence:
    @param output_file:
    @param k:
    @return:
    """
    index = process_list_to_kmer_index([query_sequence], k, docids=False, cutoff=None)
    query_hashes = [x[0] for x in index]
    data =  query_kmers(db_file, query_hashes, both_strands=True)
    with open(output_file, "w") as fh:
        for i in xrange(0, len(query_sequence)-k + 1):
            kmer = query_sequence[i:i+k]
            p0 = 0
            p1 = 0
            p2 = 0
            rkmer = get_revcomp(kmer)
            if kmer in data:
                p1 = int(data[kmer])
            if rkmer in data:
                p2 = int(data[rkmer])
            p = max(p1,p2,p0)
            fh.write("%s\t%s\t%s\n" % (kmer, i, p))


def sc_count_and_dump_kmers_for_folder(folder, output_prefix, kmers_file, k=23, mintf=None):
    """
    Count kmers, merge them, and save to tab-delimited kmers_file.
    @param folder: folder with input files
    @param output_prefix: prefix for input files
    @param kmers_file: output dump file
    @param k: kmer length
    @param mintf: minimal tf for count
    @return: None
    """
    count_kmers(input_file, output_prefix, k, mintf=mintf)
    merge_kmers(folder, input_file, input_file)
    db_file = "%s.jf" % input_file
    dump_kmers(db_file, kmers_file)

def sc_count_and_dump_kmers_for_file(fasta_file, jellyfish_data_folder, jf_db, jf_dat, k, mintf, dumpmintf):
    """
    Count kmers, merge them, and save to tab-delimited kmers_file.
    @param fasta_file: fasta file
    @param jellyfish_data_folder: output jellyfish folder
    @param jf_db: jf database output file
    @param jf_dat: jf dump output file
    @param k: kmer length
    @param mintf: minimal tf for count
    @param dumpmintf: minimal tf for dump
    @return: None
    """
    ouput_prefix = os.path.join(
            jellyfish_data_folder,
            "%s___" % jf_db,
        )
    count_kmers(fasta_file, ouput_prefix, k, mintf=mintf)
    merge_kmers(jellyfish_data_folder, ouput_prefix, jf_db)
    dump_kmers(jf_db, jf_dat, dumpmintf=dumpmintf)
    print "Sort data..."
    temp_file = jf_dat+".temp"
    data = {
        "in": jf_dat,
        "out": temp_file,
    }
    command = "sort -k2nr %(in)s > %(out)s" % data
    print command
    os.system(command)
    command = "mv %(out)s %(in)s" % data
    print command
    os.system(command)