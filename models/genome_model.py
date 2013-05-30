#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 25.12.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 


from PyExp import AbstractModel

class GenomeModel(AbstractModel):
    ''' Container for genome information.

    Dumpable attributes:

    - "genome_taxon",
    - "genome_prefix",
    - "genome_gpid",
    - "genome_acronym",
    - "genome_chromosomes",
    - "genome_contigs",
    - "genome_length",
    - "genome_mean_gc",
    - "genome_trs_all",
    - "genome_trs_3000",
    - "genome_trs_all_proc",
    - "genome_trs_3000_proc",
    - "genome_trs_all_length",
    - "genome_trs_3000_length",
    - "genome_gaps",
    - "genome_sum_gc",

    '''


    dumpable_attributes = [
        "genome_taxon",
        "genome_prefix",
        "genome_gpid",
        "genome_acronym",
        "genome_chromosomes",
        "genome_contigs",
        "genome_length",
        "genome_mean_gc",
        "genome_trs_all",
        "genome_trs_3000",
        "genome_trs_all_proc",
        "genome_trs_3000_proc",
        "genome_trs_all_length",
        "genome_trs_3000_length",
        "genome_gaps",
        "genome_sum_gc",
        ]

    def preprocess_data(self):
        ''' Preprocess data.'''
        if self.genome_trs_all_length:
            self.genome_trs_all_proc = self.genome_trs_all_length / float(self.genome_length)
        if self.genome_trs_3000_length:
            self.genome_trs_3000_proc = self.genome_trs_3000_length / float(self.genome_length)
        if not self.genome_mean_gc:
            self.genome_mean_gc = self.genome_sum_gc / self.genome_contigs