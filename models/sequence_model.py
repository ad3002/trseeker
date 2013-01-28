#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 14.06.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

from PyExp.models.abstract_model import AbstractModel
from trseeker.tools.sequence_tools import get_revcomp, get_gc, check_gapped
from trseeker.tools.sequence_tools import clear_sequence
from trseeker.tools.parsers import parse_fasta_head, parse_chromosome_name

class SequenceModel(AbstractModel):
    """ Class for sequence wrapping
    
    Dumpable attributes:
    
    - seq_gi (int)
    - seq_ref
    - seq_description
    - seq_sequence
    - seq_length (int)
    - seq_gc (float)
    - seq_revcom
    - seq_gaped (int)
    - seq_chr
    - seq_head
    - seq_start_position (int)
    - seq_end_position (int)
        
    """

    dumpable_attributes = ["seq_gi",
                           "seq_ref",
                           "seq_description",
                           "seq_sequence",
                           "seq_length",
                           "seq_gc",
                           "seq_revcom",
                           "seq_gaped",
                           "seq_chr",
                           "seq_head",
                           "seq_start_position",
                           "seq_end_position",
                           ]

    int_attributes = ["seq_gi",
                      "seq_length",
                      "seq_gaped"
                      "seq_start_position",
                      "seq_end_position", ]

    float_attributes = ["seq_gc", ]

    @property
    def length(self):
        """ Return a sequence length."""
        return self.seq_length

    @property
    def sequence(self):
        """ Return a sequence."""
        return self.seq_sequence

    @property
    def fasta(self):
        """ Return a fasta representation."""
        if not self.seq_head:
          self.seq_head = ">%s" % self.seq_ref
        return "%s\n%s\n" % (self.seq_head.strip(), self.seq_sequence.strip())

    @property
    def sa_input(self):
        """ Return a sequence left flanked by $."""
        return "%s$" % self.sequence

    @property
    def ncbi_fasta(self):
        """ Return fasta in NCBI format with gi and ref."""
        return ">gi|%s|ref|%s|%s\n%s\n" % (self.seq_gi,
                                           self.seq_ref,
                                           self.seq_description,
                                           self.seq_sequence)

    def add_sequence_revcom(self):
        """ Add reverse complement."""
        self.seq_revcom = get_revcomp(self.seq_sequence)

    def set_dna_sequence(self, title, sequence, description=None):
        """ Set sequence object."""
        self.seq_ref = title
        self.seq_description = description
        self.seq_sequence = clear_sequence(sequence)
        self.seq_length = len(self.seq_sequence)
        self.seq_gc = get_gc(self.seq_sequence)
        self.seq_gaped = check_gapped(self.seq_sequence)

    def set_ncbi_sequence(self, head, sequence):
        """ Set NCBI fasta sequence object."""
        self.seq_head = head
        self.set_dna_sequence(None, sequence)
        (self.seq_gi,
         self.seq_ref,
         self.seq_description) = parse_fasta_head(head)


        chr = parse_chromosome_name(head)
        if chr != "?":
            self.seq_chr = chr

    def set_gbff_sequence(self, head, sequence):
        """ Set NCBI fasta sequence object."""
        self.seq_head = "\t".join( [ head["gi"],
                                     head["ref"],
                                     head["description"],
                                  ])
        self.set_dna_sequence(None, sequence)
        (self.seq_gi,
         self.seq_ref,
         self.seq_description) = ( head["gi"],
                                   head["ref"],
                                   head["description"],
                                  )


        chr = parse_chromosome_name(head["description"])
        if chr != "?":
            self.seq_chr = chr
