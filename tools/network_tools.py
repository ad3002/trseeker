#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 08.09.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
Functions for network construction.

TODO: check it.
'''

from trseeker.seqio.tab_file import sc_iter_tab_file
from trseeker.models.trf_model import TRModel
import pickle
import sys
from trseeker.seqio.tr_file import read_trid2meta
import os
import array
from collections import defaultdict
try:
    import numpy
except Exception, e:
    print "Numpy import error: %s" % e
from PyExp.experiments.abstract_experiment import Timer

try:
    import graph_tool
    from graph_tool import topology
    GRAPHTOOL = True
except:
    print "WARNING: install graph_tool"
    import networkx
    GRAPHTOOL = False

def compute_network(network_file, output_file_pattern,
                    trf_large_index_file,
                    g_file, weights_file,
                    weight_to_pairs_file):
    ''' Compute network slices using graph_tool library or networkx.
    '''

    if GRAPHTOOL:
        ml_file = network_file + ".ml"
        with Timer("Create ml"):
            create_graphml(network_file, ml_file)
        with Timer("Load graph"):    
            G = load_graph(ml_file)
        with Timer("ANALYSE NETWORK"):
             analyse_network_graph_tool(G, output_file_pattern, trf_large_index_file)
    else:
        with Timer("LOAD NETWORK"):
            network_data = load_networkx(network_file)
        with Timer("INIT NETWORK"):
            G = init_graph_networkx(network_data, start=0, precise=5, 
                trf_large_index_file=trf_large_index_file)
        with Timer("ANALYSE NETWORK"):
            analyse_networkx(G, network_data, output_file_pattern, trf_large_index_file)

def load_graph(ml_file):
    ''' Load graph for graph_tools
    '''
    print "Graph loading..."
    G = graph_tool.load_graph(ml_file, file_format="xml")
    return G

def analyse_network_graph_tool(G, output_file_pattern, trf_large_index_file):
    ''' Analys graph, created by graphtool.'''
    
    print "Read meta trid"
    if trf_large_index_file:
        trid2meta = read_trid2meta(trf_large_index_file)

    print "Get edges..."
    pm = G.edge_properties["weight"]
    edges = []
    for i, e in enumerate(G.edges()):
         edges.append(e)

    print "Iterate over weights"

    last = -1
    last_component = -1
    ended = False

    np = G.vertex_properties["trid"]

    while edges:

        e = edges.pop()
        w = pm[e]
        if w == last:
            G.remove_edge(e)
            last = w
            continue
        components, hist = topology.label_components(G)
        
        N = len(hist)
        if N == last_component:
            G.remove_edge(e)
            last = w
            continue
        last_component = N

        comp_data = defaultdict(list)
        for i, c in enumerate(components.a):
            comp_data[c].append(int(np[G.vertex(i)]))

        print "Distance: %s | Edges: %s | Nodes: %s | Components: %s" % (w,
                                                                         G.num_edges(),
                                                                         i,
                                                                         N)

        if hist[0] == 1:
            ended = True
        
        output_file = output_file_pattern % (int(w), N)
        write_classification_graph_tool(output_file, comp_data, trid2meta)

        if ended:
            break

        G.remove_edge(e)
        last = w

def write_classification_graph_tool(output_file, components, trid2meta):
    ''' Save slice data to file.
    
    - components is a default dictionary: slice_id -> [trid list]
    - trid2meta is a dictionary: trid -> descripion string
    '''
    if os.path.isfile(output_file):
        os.unlink(output_file)
    with open(output_file, "w") as fw:
        for i, comp in components.items():
            for trid in comp:
                fw.write("%s\t%s" % (i, trid2meta[int(trid)]))

def write_classification_neworkx(output_file, components, trid2meta):
    ''' Save slice data to file.
    
    - components is a default dictionary: slice_id -> [trid list]
    - trid2meta is a dictionary: trid -> descripion string
    '''
    if os.path.isfile(output_file):
        os.unlink(output_file)
    with open(output_file, "w") as fw:
        for i, comp in enumerate(components):
            for trid in comp:
                fw.write("%s\t%s" % (i, trid2meta[int(trid)]))

def create_graphml(network_file, ml_file):
    ''' Create graphml xml file from tab delimited network_file.
    '''

    start = """<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns"  
    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns
     http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
<key id="d0" for="node" attr.name="trid" attr.type="double"/>
<key id="d1" for="edge" attr.name="weight" attr.type="double"/>
<graph id="G" edgedefault="undirected">
            """
    node = '<node id="n%s"><data key="d0">%s</data></node>\n'
    edge = '<edge id="e%s" source="n%s" target="n%s"><data key="d1">%s</data></edge>\n'
    end = '</graph>\n</graphml>'

    print "Load data..."
    with open(network_file, "r") as fh:
        network_data = fh.readlines()
    print "Parse data..."
    data = ""
    data += start
    seen = {}
    k = 0
    for i, line in enumerate(network_data):
        a, b, w = line.strip().split("\t")
        w = round(float(w), 4)
        if not a in seen:
            data += node % (k, a)
            seen[a] = k
            k += 1
        if not b in seen:
            data += node % (k, b)
            seen[b] = k
            k += 1
        data += edge % (i, seen[a], seen[b], w)
    data += end

    print "Save data..."
    with open(ml_file, "w") as fh:
        fh.write(data)

#
# Networkx section
#

def load_networkx(network_file):
    ''' Load network data if it was pickled previously.'''
    
    print "Network data reading..."

    a_nodes = array.array("I")
    b_nodes = array.array("I")
    weights = array.array("f")

    with open(network_file, "r") as fh:
        network_data = fh.readlines()
    print "Network data parsing..."

    for i, line in enumerate(network_data):
        a, b, w = line.strip().split("\t")
        w = float(w)
        a_nodes.append(int(a))
        b_nodes.append(int(b))
        weights.append(w)

    k = len(network_data)

    network_data = None

    print "#edges: ", k
    return (a_nodes, b_nodes, weights, k)

def init_graph_networkx(network_data, start=0, precise=1, trf_large_index_file=None):
    ''' Init graph with data.'''

    print "Graph creation..."
    G = networkx.Graph()

    print "Populate nodes..."
    for trf_obj in sc_iter_tab_file(trf_large_index_file, TRModel):
        G.add_node(int(trf_obj.trf_id))

    a_nodes, b_nodes, weight_vals, n = network_data
    print "Add edges..."
    for i in xrange(n):
        a = a_nodes[i]
        b = b_nodes[i]
        G.add_edge(a, b)
    return G

def analyse_networkx(G, network_data, output_file_pattern, trf_large_index_file):
    ''' Analize network and save slices.
    G is a networkx graph
    network_data - (a_nodes, b_nodes, weight_vals, n) parsed network data
    '''

    a_nodes, b_nodes, weight_vals, n = network_data

    print "Read meta trid"
    if trf_large_index_file:
        trid2meta = read_trid2meta(trf_large_index_file)

    print "Network analyzing..."
    last = -1
    last_component = -1
    ended = False

    k = len(weight_vals)
    while weight_vals:

        val = weight_vals.pop()
        k -= 1
        if val == last:

            G.remove_edge(a_nodes[k], b_nodes[k])
            last = val
            continue

        components = networkx.connected_components(G)
        if len(components) == last_component:
            G.remove_edge(a_nodes[k], b_nodes[k])
            last = val
            continue

        last_component = len(components)

        print "Distance: %s | Edges: %s | Nodes: %s | Components: %s" % (val,
                                                                         G.number_of_edges(),
                                                                         G.number_of_nodes(),
                                                                         len(components))

        if len(components[0]) == 1:
            ended = True
        output_file = output_file_pattern % (int(val), len(components))

        write_classification_neworkx(output_file, components, trid2meta)

        if ended:
            break

        G.remove_edge(a_nodes[k], b_nodes[k])
        last = val