#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2007-2009 Aleksey Komissarov ( ad3002@gmail.com )
# All rights reserved.
#
# This software is licensed as described in the file COPYING, which
# you should have received as part of this distribution.
"""
Python BLAST wrapper.

- blastn(database, input, output)
- create_db(fasta_file, verbose=False, base_name=False, title=False)
- bl2seq(input1, input2, output)

Used settings

>>> location = settings["blast_settings"]["blast_location"]
>>> b = settings["blast_settings"]["blast_b"]
>>> e = settings["blast_settings"]["blast_e"]
>>> v = settings["blast_settings"]["blast_v"]
>>> OS = settings["trseeker"]["os"]

"""


from trseeker.seqio.tab_file import sc_iter_tab_file
from trseeker.models.blast_model import BlastResultModel
from trseeker.models.trf_model import TRModel
import os
import re

from trseeker.settings import load_settings
from trseeker.seqio.tr_file import get_all_trf_objs

settings = load_settings()
location = settings["blast_settings"]["blast_location"]
if not location:
    location = ""
b = settings["blast_settings"]["blast_b"]
e = settings["blast_settings"]["blast_e"]
v = settings["blast_settings"]["blast_v"]
OS = settings["trseeker"]["os"]

def blastn(database, query, output, e_value=None):
    """ Blastn. 

    Input format parametr:

    7 qseqid qgi qacc sseqid means sgi sacc qstart qend sstart send evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe

    """

    if OS == "win":
        progr_name = "blastn.exe"
    else:
        progr_name = "blastn"

    # TODO: parameters
    format = '"7 qseqid qgi qacc sseqid sallseqid means sgi sallgi sacc sallacc qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe"'

    format = '"7 qseqid qgi qacc sseqid means sgi sacc qstart qend sstart send evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe"'

    if not e_value:
        e_value = e

    string = '%s%s -query %s -task blastn -db %s -out %s -evalue %s -word_size 10 -outfmt %s -dust no -soft_masking "false" -max_target_seqs %s -num_descriptions %s -num_threads 8 -num_alignments %s' % (location,
                          progr_name,
                          query,
                          database,
                          output,
                          e_value,
                          format,
                          b,
                          v,
                          b,
                          )
    print string
    os.system(string)

    # Examples:
    #    blastn -query HTT_gene -task megablast -db hs_chr -db_soft_mask 30 -outfmt 7 -out HTT_megablast_mask.out -num_threads 4
    #    $ echo 1786181 | ./blastn -db ecoli -outfmt "7 qacc sacc evalue qstart qend sstart send"


def create_db(fasta_file, output, verbose=False, title=None):
    """ Create BLAST database.

    >>> # string
    >>> "%s%s -in %s -hash_index -dbtype nucl -parse_seqids -out %s" % (location, progr_name, fasta_file, output)

    """

    # TODO: parameters

    if verbose:
        print "Format fasta file for BLAST search ...";

    if OS == "win":
        progr_name = "makeblastdb.exe"
    else:
        progr_name = "makeblastdb"

    string = "%s%s -in %s -hash_index -dbtype nucl -parse_seqids -out %s" % (location, progr_name, fasta_file, output)
    if title:
        string += " -title %s" % title
    if verbose:
        print string
    os.system(string)
    if verbose:
        print "Format fasta file for BLAST search complete!"

    # formatdb.exe -i E:\home\ad3002\work\mouse_wgs\fa\name.fa -p F -o T -V T -n mouse_wgs -t mouse_wgs
    # makeblastdb -in hs_chr.fa    hs_chr -title "Human chromosomes, Ref B37.1"

def alias_tool(dblist, output, title):
    """ Create alias database """

    if OS == "win":
        progr_name = "blastdb_aliastool.exe"
    else:
        progr_name = "blastdb_aliastool"
    string = "%s%s -dblist %s -dbtype nucl -out %s -title %s" % (location, progr_name, dblist, output, title)
    print string
    os.system(string)

    # blastdb_aliastool -dblist "nematode_mrna nematode_genomic" -dbtype nucl -out nematode_all -title "Nematode RefSeq mRNA + Genomic"

def bl2seq(input1, input2, output):
    """ Blast two seq."""

    if OS == "win":
        _location = settings["blast_settings"]["blast_location_WIN"]
    else:
        _location = settings["blast_settings"]["blast_location_NIX"]

    # TODO: parameters
    string = "perl %s -p blastn -i %s -j %s -o %s -F F -W 5 -D 1 -e %s" % (
                                          _location,
                                          input1,
                                          input2, 
                                          output, 
                                          e)
    print _location
    print string[:100]
    os.system(string)
    exit()

def create_db_for_genome(file_pattern=None,
                         chromosome_list=None,
                         output=None,
                         title=None
                          ):
    ''' Create BLAST DB for genome.

    Example:

    >>> create_db_for_genome(file_pattern="E:/eukaryota_genomes/mus_musculus/fasta/%s.fsa_nt",
    ...                 chromosome_list=["wgs.AAHY.1", "wgs.AAHY.2", "wgs.AAHY.3", "wgs.AAHY.4", "wgs.AAHY.5", "wgs.AAHY.6", "wgs.CAAA"],
    ...                 output="mus_musculus_wgsa",
    ...                 title="Mouse WGSAs"
    ...                   )

    '''
    files = []
    for x in chromosome_list:

        fasta_file = file_pattern % x
        print "DB construction: ", fasta_file

        create_db(fasta_file, fasta_file, verbose=True, title=fasta_file.split("/")[-1])

        files.append(fasta_file)
    dblist = '"' + " ".join(files) + '"'
    title = '"%s"' % title
    alias_tool(dblist, output, title)


def get_gi_list(gi_score_file, score_limit=90):
    ''' Function return GI list from given score file. '''
    result = []
    with open(gi_score_file, "rb") as fh:
        for line in fh:
            data = line.strip().split("\t")
            if float(data[1]) > score_limit:
                result.append(data[0])
    return result

def get_all_gi_from_blast(blast_file, mode="gi"):
    ''' Get gi or ref -> hits list from blast output.'''

    result = {}
    for blast_obj in sc_iter_tab_file(blast_file, BlastResultModel, skip_starts_with="#"):

        if mode == "gi":
            id = blast_obj.subject_gi
        elif mode == "ref":
            id = blast_obj.subject_ref
        else:
            raise Exception

        result.setdefault(id, 0)
        result[id] += 1

    result = [ (id, result[id]) for id in result.keys() ]
    result.sort(reverse=True, key=lambda x: x[1])
    return result

def get_all_blast_obj_from_blast(blast_file, mode="ref"):
    ''' Get all gi to blast_obj list dictionary from blast output.'''
    result = {}
    for blast_obj in sc_iter_tab_file(blast_file, BlastResultModel, remove_starts_with="#"):
        if mode == "gi":
            id = blast_obj.subject_gi
        elif mode == "ref":
            id = blast_obj.subject_ref
        else:
            raise Exception
        result.setdefault(id, [])
        result[id].append(blast_obj)
    return result
        
def bl2seq_search_for_trs(trf_large_file, annotation_bl2seq_folder, temp_file):
    ''' Pairwise search of TRs with bl2seq.
    '''
    # read trf_objs
    print "Read TRs..."
    trf_objs = get_all_trf_objs(trf_large_file)
    N = len(trf_objs)
    for i in xrange(N):
        for j in xrange(i, N):
            print "Searching:", i, j, "from", N
            output_file = "%s.%s.blast" % (trf_objs[i].trf_id, trf_objs[j].trf_id)
            blast_output_file = os.path.join(annotation_bl2seq_folder, output_file)
            input1 = trf_objs[i].trf_array
            input2 = trf_objs[j].trf_array
            input_file1 = temp_file + ".1.fasta"
            input_file2 = temp_file + ".2.fasta"
            with open(input_file1, "w") as fh:
                fh.write(">%s\n%s" % (trf_objs[i].trf_id, input1))
            with open(input_file2, "w") as fh:
                fh.write(">%s\n%s" % (trf_objs[j].trf_id, input2))
            bl2seq(input_file1, input_file2, blast_output_file)

def blastn_search_for_trs(trf_large_file, db, annotation_self_folder, temp_file, skip_by_family=None, is_huge_alpha=False):
    ''' Search TRs in given DB.
    '''

    # count files
    for n, u in enumerate(sc_iter_tab_file(trf_large_file, TRModel)):
        pass
    if is_huge_alpha:
        alpha_sets = {}
    for i, u in enumerate(sc_iter_tab_file(trf_large_file, TRModel)):
        print "%s/%s" % (i, n)
        
        if skip_by_family:
            if u.trf_family in skip_by_family:
                continue
        if u.trf_family == "ALPHA":
            continue
        
        blast_output_file = os.path.join(annotation_self_folder, "%s.blast" % u.trf_id)
        if os.path.isfile(blast_output_file):
            with open(blast_output_file) as fh:
                    data = fh.read()
            if data.startswith("ALPHA"):
                continue
            
        seen = False
        if is_huge_alpha:
            for key, tr_set in alpha_sets.items():
                if u.trf_id in tr_set:
                    # print "SEEN ALPHA", key, u.trf_id
                    with open(blast_output_file, "w") as fh:
                        fh.write("ALPHA\t%s" % key)
                    seen = True
        if seen:
            continue

        
        if os.path.isfile(blast_output_file) and os.stat(blast_output_file).st_size != 0:
            trids = _get_trids_from_blast_file(blast_output_file)
            alpha_sets[u.trf_id] = trids
            print "ADDED ALPHA by %s with length %s" % (u.trf_id, len(trids)) 
            continue

        with open(temp_file, "w") as fh:
            fh.write(">%s\n%s" % (u.trf_id, u.trf_array))
        blastn(db, temp_file, blast_output_file)
        if os.path.isfile(temp_file):
            os.remove(temp_file)

        if is_huge_alpha:
            trids = _get_trids_from_blast_file(blast_output_file)
            alpha_sets[u.trf_id] = trids
            print "ADDED ALPHA by %s with length %s" % (u.trf_id, len(trids)) 
            
def _get_trids_from_blast_file(blast_output_file):
    '''
    '''
    result = set()
    gi2obj = get_all_blast_obj_from_blast(blast_output_file)
    for gi, hits in gi2obj.items():
        for hit in hits:
            if hit.score >= 90:
                result.add(int(hit.subject_ref))
    return list(result)
            