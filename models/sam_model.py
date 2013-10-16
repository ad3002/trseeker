#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 28.11.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

from PyExp import AbstractModel

class SAMModel(AbstractModel):
    """ Class for SAM data wrapping.

    Dumpable attributes:

    - QNAME
    - FLAG
    - RNAME
    - POS
    - MAPQ
    - CIGAR
    - RNEXT
    - PNEXT
    - TLEN
    - SEQ
    - QUAL
    - features
    
    """

    dumpable_attributes = [
        "QNAME",
        "FLAG",
        "RNAME",
        "POS",
        "MAPQ",
        "CIGAR",
        "RNEXT",
        "PNEXT",
        "TLEN",
        "SEQ",
        "QUAL",
        "features",
    ]
    
    int_attributes = [
        "FLAG",
        "POS",
        "MAPQ",
        "PNEXT",
        "TLEN",
    ]

    @property
    def get_fragment_length(self):
        return self.POS - self.PNEXT