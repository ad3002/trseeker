#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2007-2009 Aleksey Komissarov ( ad3002@gmail.com )
# All rights reserved.
#
# This software is licensed as described in the file COPYING, which
# you should have received as part of this distribution.
"""
BLAST wrapper.
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
    """ Blast fasta file versus blast database.
    Available output format parameters:
        7 qseqid qgi qacc sseqid means sgi sacc qstart qend sstart send evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe
    Command examples:
        blastn -query HTT_gene -task megablast -db hs_chr -db_soft_mask 30 -outfmt 7 -out HTT_megablast_mask.out -num_threads 4
        $ echo 1786181 | ./blastn -db ecoli -outfmt "7 qacc sacc evalue qstart qend sstart send"
    @todo: do something with parameters
    @param database: database created with makeblastdb
    @param query: fasta file
    @param output: output file
    @param e_value: e value, e.g. 1e-20
    @return: None
    """
    if OS == "win":
        program_name = "blastn.exe"
    else:
        program_name = "blastn"
    format_string = '"7 qseqid qgi qacc sseqid means sgi sacc qstart qend sstart send evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe"'

    if not e_value:
        e_value = e

    data = {
        'location': location,
        'program_name': program_name,
        'query': query,
        'database': database,
        'output': output,
        'e_value': e_value,
        'format_string': format_string,
        'b': b,
        'v': v,
    }

    string = '%(location)s%(program_name)s -query %(query)s -task blastn -db %(database)s ' +\
             '-out %(output)s -evalue %(e_value)s -word_size 10 -outfmt %(format_string)s ' +\
             '-dust no -soft_masking "false" -max_target_seqs %(b)s -num_descriptions %(v)s ' +\
             '-num_threads 8 -num_alignments %(b)s' % data

    print string
    os.system(string)


def create_db(fasta_file, output, verbose=False, title=None):
    """  Create BLAST database.
    Command example:
        "%s%s -in %s -hash_index -dbtype nucl -parse_seqids -out %s" % (location, progr_name, fasta_file, output)
        formatdb.exe -i E:\home\ad3002\work\mouse_wgs\fa\name.fa -p F -o T -V T -n mouse_wgs -t mouse_wgs
        makeblastdb -in hs_chr.fa    hs_chr -title "Human chromosomes, Ref B37.1"
    @todo: do something with parameters
    @param fasta_file: fasta file
    @param output: output database path
    @param verbose: verbose
    @param title: database title
    @return: None
    """
    if verbose:
        print "Format fasta file for BLAST search ..."

    if OS == "win":
        program_name = "makeblastdb.exe"
    else:
        program_name = "makeblastdb"

    string = "%s%s -in %s -hash_index -dbtype nucl -parse_seqids -out %s" % (location, program_name, fasta_file, output)
    if title:
        string += " -title %s" % title
    if verbose:
        print string
    os.system(string)
    if verbose:
        print "Format fasta file for BLAST search complete!"


def alias_tool(dblist, output, title):
    """ Create alias database
    Command example:
        blastdb_aliastool -dblist "nematode_mrna nematode_genomic" -dbtype nucl -out nematode_all -title "Nematode RefSeq mRNA + Genomic"
    @param dblist: list of database paths
    @param output: output database path
    @param title: alias database title
    @return: None
    """
    if OS == "win":
        program_name = "blastdb_aliastool.exe"
    else:
        program_name = "blastdb_aliastool"
    string = "%s%s -dblist %s -dbtype nucl -out %s -title %s" % (location, program_name, dblist, output, title)
    print string
    os.system(string)


def bl2seq(input1, input2, output):
    """
    Compare two sequences with blast.
    @todo: do something with parameters
    @param input1: file with fasta
    @param input2: file with fasta
    @param output: output file
    @return: None
    """
    if OS == "win":
        _location = settings["blast_settings"]["blast_location_WIN"]
    else:
        _location = settings["blast_settings"]["blast_location_NIX"]

    string = "perl %s -p blastn -i %s -j %s -o %s -F F -W 5 -D 1 -e %s" % (
                                                        _location,
                                                        input1,
                                                        input2,
                                                        output,
                                                        e)
    print string
    os.system(string)


def create_db_for_genome(file_pattern=None, chromosome_list=None, output=None, title=None):
    """
    Create BLAST database for genome.
    Example:
        create_db_for_genome(file_pattern="E:/eukaryota_genomes/mus_musculus/fasta/%s.fsa_nt",
                     chromosome_list=["wgs.AAHY.1", "wgs.AAHY.2", "wgs.AAHY.3", "wgs.AAHY.4", "wgs.AAHY.5", "wgs.AAHY.6", "wgs.CAAA"],
                     output="mus_musculus_wgsa",
                     title="Mouse WGSAs")
    @param file_pattern: path pattern
    @param chromosome_list: list of chromosomes for file pattern
    @param output: final database path
    @param title: final database title
    @return: None
    """
    files = []
    for x in chromosome_list:
        fasta_file = file_pattern % x
        print "Database construction: ", fasta_file
        create_db(fasta_file, fasta_file, verbose=True, title=fasta_file.split("/")[-1])
        files.append(fasta_file)
    dblist = '"' + " ".join(files) + '"'
    title = '"%s"' % title
    alias_tool(dblist, output, title)


def get_gi_list(gi_score_file, score_limit=90):
    """
    Get GI list from given score file.
    @param gi_score_file: GI score file
    @param score_limit: get only GI with score greater than score limit (default 90)
    @return: list of GI
    """
    result = []
    with open(gi_score_file, "rb") as fh:
        for line in fh:
            data = line.strip().split("\t")
            if float(data[1]) > score_limit:
                result.append(data[0])
    return result


def get_all_gi_from_blast(blast_file, mode="gi"):
    """
    Get all GI from blast results file.
    @param blast_file: blast file
    @param mode: parsing mode GI or ref based
    @return: a dictionary of gi or ref to number of hist
    """
    result = {}
    for blast_obj in sc_iter_tab_file(blast_file, BlastResultModel, remove_starts_with="#"):
        if mode == "gi":
            _id = blast_obj.subject_gi
        elif mode == "ref":
            _id = blast_obj.subject_ref
        else:
            raise Exception("Unknown mode")
        result.setdefault(_id, 0)
        result[_id] += 1

    result = [(_id, result[_id]) for _id in result.keys()]
    result.sort(reverse=True, key=lambda x: x[1])
    return result


def get_all_blast_obj_from_blast(blast_file, mode="ref"):
    """
    Get all GI to blast_obj list dictionary from blast output.
    @param blast_file:
    @param mode:
    @return:
    """
    result = {}
    for blast_obj in sc_iter_tab_file(blast_file, BlastResultModel, remove_starts_with="#"):
        if mode == "gi":
            _id = blast_obj.subject_gi
        elif mode == "ref":
            _id = blast_obj.subject_ref
        else:
            raise Exception
        result.setdefault(_id, [])
        result[id].append(blast_obj)
    return result


def bl2seq_search_for_trs(trf_large_file, annotation_bl2seq_folder, temp_file):
    """
    Pairwise search of TRs with bl2seq.
    @param trf_large_file:
    @param annotation_bl2seq_folder:
    @param temp_file:
    @return:
    """
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
    """
    Search TRs in given DB.
    @param trf_large_file:
    @param db:
    @param annotation_self_folder:
    @param temp_file:
    @param skip_by_family:
    @param is_huge_alpha:
    @return:
    """
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
    """
    @param blast_output_file:
    @return:
    """
    result = set()
    gi2obj = get_all_blast_obj_from_blast(blast_output_file)
    for gi, hits in gi2obj.items():
        for hit in hits:
            if hit.score >= 90:
                result.add(int(hit.subject_ref))
    return list(result)