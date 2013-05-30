#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 26.12.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
Functions related to Repbase datasets.
'''
import os
from trseeker.tools.sequence_tools import clear_sequence
from PyExp import sc_iter_filepath_folder
from trseeker.seqio.fasta_file import sc_iter_fasta

def _join_repbase_fastas(fasta_file, fasta_output, index_file_name=None, file_name=None, start_id=None, increment=False, repbase_names_fix=False):
    ''' Function join all Repbase fasta files in one huge fasta; reformat headers for compatibility with NCBI tools
    '''
    if not increment:
        if os.path.isfile(fasta_output):
            os.unlink(fasta_output)

    with open(fasta_output, "a") as fa_fw:
        if start_id:
            i = start_id
        else:
            i = 0
        for seq_obj in sc_iter_fasta(fasta_file):
            i += 1
            head = seq_obj.seq_head
            if file_name and repbase_names_fix:
                file_name = file_name.split(".")[0]
                head_info = "%s:%s" % (file_name, head.replace("  ", " ").replace(" ", "_").replace(">", "").replace("(", "").replace(")", "").replace("\t", ":").strip())
            else:
                head_info = head.strip().replace(">", "")

            fa_fw.write(">%s\n%s\n" % (head_info, seq_obj.sequence))
    return i

def join_repbase_files(input_folder, output_file):
    ''' Function join all Repbase fasta files in one huge fasta; reformat headers for compatibility with NCBI tools.'''
    start_id = 0
    for input_file in sc_iter_filepath_folder(input_folder, mask="."):
        if not input_file.endswith("ref"):
            continue
        start_id = _join_repbase_fastas(input_file,
                                   output_file,
                                   file_name=os.path.split(input_file)[-1],
                                   start_id=start_id,
                                   increment=True,
                                   repbase_names_fix=True)
