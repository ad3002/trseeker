#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 05.06.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

"""
Classes:

- TRFFileIO(AbstractBlockFileIO)

"""
import re, os
from trseeker.models.sam_model import SAMModel
from trseeker.seqio.tab_file import TabDelimitedFileIO
from PyExp import AbstractFileIO
from PyExp import WizeOpener

class SAMFileIO(TabDelimitedFileIO):
    """ 
    """

    def __init__(self, *args, **kwargs):
        """ 
        """
        super(TabDelimitedFileIO, self).__init__(*args, **kwargs)
        self.headers = {}

    def read_online(self, file_name):
        """ Overrided. Yield items online from data from input_file."""
        with open(input_file) as fh:
            for i, line in enumerate(fh):
                if line.startswith("@"):
                    self.headers.append()
                else:
                    break
        fields = SAMModel().dumpable_attributes
        for data in csv.DictReader(fh, fieldnames=fields, delimiter='\t', quoting=csv.QUOTE_NONE):
            obj = data_type()
            obj.set_with_dict(data)
            yield obj

def sc_sam_reader(sam_file):
    """
    """
    reader = SAMFileIO()
    