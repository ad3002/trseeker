#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 07.09.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
Functions related to assembler statistics
'''

def get_N50(lengths):
    '''
    Get (N50 contig length, N50) statistics for list of contigs lengths
    '''
    lengths.sort(reverse=True)
    total = sum(lengths)
    n50 = 0
    l50 = 0
    for x in lengths:
        l50 += 1
        n50 += x
        if n50 >= total/2:
            return x, l50