#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 05.06.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
SRA files.
TODO: implement this.
TODO: separate model from reader
'''

from trseeker.tools.sequence_tools import get_revcomp, get_gc
from trseeker.tools.sequence_tools import clear_sequence
from Bio import pairwise2

class FastqObj(object):

    def __init__(self, head, seq, strain, qual, phred33=False):
        self.head = head.strip()
        self.seq = seq.lower().strip()
        self.qual = qual.strip()
        self.strain = strain.strip()
        self.id = head.strip().split()[0].replace("@","")
        # trimmeing input
        self.trimmed_seq = self.seq.strip()
        self.trimmed_qual = qual.strip()
        if phred33:
            self.qv = [ord(x)-33 for x in qual.strip()]
        # trimming results
        self.adaptor_positions = []
        self.adapter_contamination = None
        self.qual_junk = None
        self.status = True
        self.parts = []
        

    @property
    def fastq(self):
        return "%s\n%s\n%s\n%s\n" % (
                    self.head,
                    self.seq,
                    self.strain,
                    self.qual,
            )

    def fastq_with_error(self, error):
        return "@%s__%s\n%s\n%s\n%s\n" % (
                    error,
                    self.head[1:],
                    self.seq,
                    self.strain,
                    self.qual,
            )

    @property
    def trimmed_fastq(self):
        return "%s%s\n%s%s\n" % (
                    self.head,
                    self.trimmed_seq,
                    self.strain,
                    self.trimmed_qual,
            )

    @property
    def sequence(self):
        return self.seq.strip().lower()

    @property
    def fasta(self):
        return ">%s\n%s\n" % (
                    self.head.strip(),
                    self.seq.strip(),
                    )

    @property
    def gc(self):
        return get_gc(self.seq)

    @property
    def length(self):
        return len(self.seq)

    def trim(self):
        '''
        '''
        raise NotImplemented
                    
    def trim_by_quality_cutoff(self, cutoff=30):
        '''
        '''
        for i, q in enumerate(self.qv):
            if q < cutoff:
                self.trimmed_seq = self.trimmed_seq[:i]
                self.trimmed_qual = self.trimmed_qual[:i]
                self.qual_junk = self.trimmed_seq[i:]
                break

    def trim_exact_adaptor(self, adaptor, verbose=False):
        '''
        '''
        raise NotImplemented
       
    def trim_inexact_adaptor(self, adaptor, verbose=False, num_errors=2):
        '''
        '''
        raise NotImplemented

class PERun(object):
    pass

def fastq_reader(fastq_file, phred33=False):
    with open(fastq_file) as fh:
        while True:
            try:
                head = fh.readline()
                seq = fh.readline()
                strain = fh.readline()
                qual = fh.readline()
                fastq_obj = FastqObj(head, seq, strain, qual, phred33=phred33)
                yield fastq_obj
            except:
                break