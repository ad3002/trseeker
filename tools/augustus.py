#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2007-2022 Aleksey Komissarov ( ad3002@gmail.com )
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

import os, shutil, tempfile
from multiprocessing import Pool
from trseeker.seqio.fasta_file import sc_iter_fasta


def augustus_search_by_splitting(fasta_file, output_file, threads=30, wdir="."):
    """ Augustus search by splitting on fasta file in files.
    """
    folder_path = tempfile.mkdtemp(dir=wdir)
    
    ### 1. Split chromosomes into temp file
    total_length = 0
    next_file = 0
    for i, seq_obj in enumerate(sc_iter_fasta(fasta_file)):
        print(seq_obj.header)
        file_path = os.path.join(folder_path, "%s.fa" % next_file)
        with open(file_path, "a") as fw:
            fw.write(">%s\n%s\n" % (seq_obj.header, seq_obj.sequence))
        total_length += len(seq_obj.sequence)
        if total_length > 100000:
            next_file += 1
            total_length = 0

    ### 2. Run Augustus

    current_dir = os.getcwd()

    os.chdir(folder_path)
    
    command = "ls %s | grep fa | xargs -P %s -I {} sh -c 'augustus --species=human --strand=both --genemodel=partial --gff3=on --progress=true {} --uniqueGeneId=true > {}.augustas'" % (folder_path, threads)
    print(command)
    os.system(command)

    
    ### 3. Aggregate data

    command = "cat %s/*.augustus > %s" % (folder_path, output_file)
    print(command)
    os.system(command)

    os.chdir(current_dir)

    ### 4. Remove temp folder
    input("Remove: %s ?" % folder_path)
    shutil.rmtree(folder_path)