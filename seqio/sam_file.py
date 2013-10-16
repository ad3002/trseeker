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

    def read(self, file_name)   

def sc_sam_reader(sam_file):
    """
    """
    reader = SAMFileIO()
    