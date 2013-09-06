#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 05.06.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
"""
Classes:

- FastaFileIO(AbstractBlockFileIO)
  
"""
from trseeker.seqio.block_file import AbstractBlockFileIO
from trseeker.models.sequence_model import SequenceModel

class FastaFileIO(AbstractBlockFileIO):
    """  Working with multi fasta files, 
    where each block starts with '>' token.
    
    Overrided public methods:
    
    - __init__(self)
    
    Inherited public properties:
    
    - data  - iterable data, each item is tuple (head, body)
    - N     - a number of items in data
    
    Inherited public methods:
    
    - read_from_file(self, input_file)
    - read_online(self, input_file) ~> item
    - get_block_sequence(self, head_start, next_head, fh)
    - get_blocks(self, token, fh)
    - gen_block_sequences(self, token, fh)
    - read_from_db(self, db_cursor)
    - write_to_file(self, output_file)
    - write_to_db(self, db_cursor)
    - read_as_iter(self, source)
    - iterate(self) ~> item of data 
    - do(self, cf, args) -> result
    - process(self, cf, args)
    - clear(self)
    - do_with_iter(self, cf, args) -> [result,]
    - process_with_iter(self, cf, args)
        
    """

    def __init__(self):
        """ Overrided. Hardcoded start token."""
        token = ">"
        super(FastaFileIO, self).__init__(token)

def sc_iter_fasta(file_name):
    """ Iter over fasta file."""

    reader = FastaFileIO()
    for (head, sequence, start, next) in reader.read_online(file_name):
        seq_obj = SequenceModel()
        seq_obj.set_ncbi_sequence(head, sequence)
        yield seq_obj

def fasta_reader(file_name):
    """ Synonym for  sc_iter_fasta.
    """
    return  sc_iter_fasta(file_name)

def sc_iter_fasta_simple(file_name):
    """ Iter over fasta file."""

    reader = FastaFileIO()
    for (head, sequence, start, next) in reader.read_online(file_name):
        seq_obj = SequenceModel()
        seq_obj.set_ncbi_sequence(head, sequence)
        yield seq_obj.seq_gi, seq_obj.sequence

def save_all_seq_with_exact_substring(fasta_file, substring, output_file):
    ''' Save all fasta sequences with exact match with given substring in output_file.
    '''
    with open(output_file, "w") as fh:
        for seq_obj in sc_iter_fasta(fasta_file):
            if substring in seq_obj.seq_sequence:
                fh.write(seq_obj.fasta)

def sort_fasta_file_by_length(file_name):
    ''' Sort by length and save fasta file
    '''
    objs = []
    for i, seq_obj in enumerate(sc_iter_fasta(fasta_file)):
        objs.append(seq_obj)
    objs.sort(key=lambda x: x.length)
    with open(fasta_file, "w") as fh:
        for obj in objs:
            fh.write(obj.fasta)
