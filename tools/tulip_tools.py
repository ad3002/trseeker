#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 07.09.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
Functions related to tulip graph format.
'''

import csv

def format_distances_to_tulip(distance_file, tulip_file_output, cutoff=90):
    """ Format and write network data to tulip format from distance_file with given cutoff (defaul=90).
    
    - distance_file: (node_a, node_b, distance)
    - tulip_file_output: network in tulip format
    - cutoff: distance cutoff
    """

    NODES = set()
    EDGES = set()

    with open(distance_file) as fh:
        for node_a, node_b, distance in csv.reader(fh, delimiter='\t'):
            NODES.add(int(node_a))
            NODES.add(int(node_b))
            EDGES.add((node_a, node_b, float(distance)))

    with open(tulip_file_output, "w") as fw:

        # write header
        fw.write("(tlp \"2.0\";\n")
        # write nodes
        node_str = "(nodes"
        for node in NODES:
            node_str = node_str + " " + str(node)
        node_str = node_str + ")\n"
        fw.write(node_str)
        # write edges
        i = 0
        for edge in EDGES:
            if edge[2] > cutoff:
                i += 1
                edge_str = "(edge %s %s %s)\n" % (i, edge[0], edge[1])
                fw.write(edge_str)
        fw.write(")")
