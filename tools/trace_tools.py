#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 28.11.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com
'''
Functions related to Trace data.
'''

import os
from trseeker.seqio.tab_file import sc_iter_simple_tab_file
from trseeker.seqio.fasta_file import sc_iter_fasta

def unclip_trace_file(fasta_file, clip_file, uncliped_file):
    ''' Unclip. Remove vector flanks from Trace data.'''

    id2clip = {}


    if os.path.isfile(clip_file):
        for (id, start, end) in sc_iter_simple_tab_file(clip_file):
            if id == "TI":
                continue
            id2clip[int(id)] = (int(start), int(end))

    result = []
    for seq_obj in sc_iter_fasta(fasta_file):
        id = int(seq_obj.seq_gi)
        seq = seq_obj.sequence
        if id in id2clip:
            seq = seq_obj.sequence[id2clip[id][0]: id2clip[id][1]]
        result.append(">%s\n%s\n" % (id, seq))

    with open(uncliped_file, "w") as fh:
        fh.writelines(result)
