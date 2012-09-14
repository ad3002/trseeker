#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 01.05.2012
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

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
import os
import re
from trseeker.seqio.block_file import AbstractBlockFileIO
from trseeker.models.sequence_model import SequenceModel
from trseeker.tools.sequence_tools import clear_sequence
from PyExp.readers.abstract_reader import sc_iter_filename_folder

class GbffFileIO(AbstractBlockFileIO):
    """ Working with multi fasta files, where each block starts with '>' token.
        
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
        token = "LOCUS       "
        super(GbffFileIO, self).__init__(token)

def sc_iter_gbff(file_name):
    """ Iter over gbff file."""

    reader = GbffFileIO()
    for (head, gb_body, head_start, next_head) in reader.read_online(file_name):

        head_data = {}
        head = head.strip()

        try:
            sequence = re.findall("ORIGIN(.*?)//", gb_body, re.S)[0]
            sequence = clear_sequence(sequence)
        except:
            sequence = ""
        position = (head_start, next_head)
        head_data["ref"] = head.split("       ")[1]
        try:
            head_data["gi"] = re.findall("GI:(\d*)\s", gb_body, re.S)[0]
        except:
            head_data["gi"] = head_data["ref"]

        try:
            head_data["description"] = re.findall("DEFINITION(.*?)ACCESSION", gb_body, re.S)[0]
            head_data["description"] = re.sub("\s+", " ", head_data["description"])
            head_data["description"] = head_data["description"].strip()
        except:
            head_data["description"] = "unknown"
        # Initiate sequence object
        seq_obj = SequenceModel()
        seq_obj.set_gbff_sequence(head_data, sequence)
        yield seq_obj

def sc_iter_gbff_simple(file_name):
    """ Iter over gbff file."""

    reader = GbffFileIO()
    for seq_obj in  sc_iter_gbff(file_name):
        yield seq_obj.seq_gi, seq_obj.sequence


def sc_parse_gbff_in_folder(gbbf_folder, fasta_folder, fasta_postfix, mask):
    ''' Parse all gbff files in given gbbf_folder according to mask, write result
    to fasta_folder with corresponding fasta_postfix.
    '''
    for file_name in sc_iter_filename_folder(gbbf_folder, mask):

    	print "Convert file: ", file_name

    	fasta_name = ".".join( file_name.split(".")[:-1] )

        file_name = os.path.join(gbbf_folder, file_name)
        


        fasta_output = os.path.join(fasta_folder, 
                                    "%s.%s" % (fasta_name,fasta_postfix)
                                    )
        print file_name, fasta_output
        if os.path.isfile(fasta_output):
            os.remove(fasta_output)
        with open(fasta_output, "ab") as fw:
            for seq_obj in sc_iter_gbff(file_name):
                fw.write(seq_obj.ncbi_fasta)
