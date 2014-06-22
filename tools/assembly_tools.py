#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 07.09.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
"""
Functions related to assembly statistics
"""


def get_n50(lengths):
    """ Get (N50 contig length, N50, shortest contig, longest contig) statistics for list of contigs lengths.
    @param lengths: a list of contig lengths
    @return (N50, L50, shortest contig, longest contig)
    """
    if not lengths:
        return 0, 0, 0, 0
    lengths.sort(reverse=True)
    total = sum(lengths)
    n50 = 0
    l50 = 0
    shortest_seq = min(lengths)
    longest_seq = max(lengths)
    for x in lengths:
        l50 += 1
        n50 += x
        if n50 >= total/2:
            return x, l50, shortest_seq, longest_seq