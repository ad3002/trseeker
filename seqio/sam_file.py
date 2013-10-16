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
import csv


class SAMFileIO(TabDelimitedFileIO):
    """ 
    """

    def __init__(self, *args, **kwargs):
        """ 
        """
        super(TabDelimitedFileIO, self).__init__(*args, **kwargs)
        self.headers = []

    def read_online(self, file_name):
        """ Overrided. Yield items online from data from input_file."""

        def skip_comments(iterable):
            for line in iterable:
                if not line.startswith('@'):
                    yield line

        with open(file_name) as fh:
            for i, line in enumerate(fh):
                if line.startswith("@"):
                    self.headers.append(line)
                else:
                    break

        fields = SAMModel().dumpable_attributes

        with open(file_name) as fh:
            for data in csv.DictReader(skip_comments(fh), fieldnames=fields, delimiter='\t', quoting=csv.QUOTE_NONE):
                k, t, v = data["features"].split(":")
                if t == "i":
                    t = int
                _features = {
                    k: t(v),
                }
                if data.has_key(None):
                    for line in data[None]:
                        k, t, v = line.split(":")
                        if t == "i":
                            t = int
                        else:
                            t = str
                        _features[k] = t(v)
                    del data[None]
                data["features"] = _features
                obj = SAMModel()
                obj.set_with_dict(data)
                yield obj

def sc_sam_reader(sam_file):
    """
    """
    reader = SAMFileIO()
    for sam_obj in reader.read_online(sam_file):
        yield sam_obj