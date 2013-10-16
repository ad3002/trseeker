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
from PyExp import AbstractFileIO

class SAMFileIO(AbstractFileIO):
    """          
    """
    def __init__(self):
        """ 
        """
        super(AbstractFileIO, self).__init__()

    def get_opener(self):
        return WizeOpener

    def read_from_file(self, input_file):
        """ Read data from given input_file."""
        with WizeOpener(input_file) as fh:
            self._data = fh.readlines()

    def read_online(self, input_file):
        """ Yield items from data online from input_file."""
        with WizeOpener(input_file) as fh:
            for item in fh:
                yield item

    def read_from_db(self, db_cursor):
        """ Read data from database cursor."""
        for item in db_cursor:
            yield item

    def read_from_mongodb(self, table, query):
        """ Read data online from mongodb."""
        cursor = table.find(query)
        n = cursor.count(query)
        start = 0
        limit = 2000
        end = start + limit
        while True:
            for x in cursor[start:end]:
                yield x
            start = end
            end = end + limit
            if start > n:
                break

    def update_mongodb(self, table, what, wherewith):
        table.update(what, wherewith, False, True)

    def write_to_file(self, output_file):
        """ Write data to given output_file."""
        with WizeOpener(output_file, "w") as fh:
            fh.writelines(self._data)

    def write_to_db(self, db_cursor):
        """ Write data with given database cursor."""
        raise NotImplementedError

    def write_to_mongodb(self, table, item):
        table.insert(item)

    def read_as_iter(self, source):
        """ Read data from iterable source."""
        for item in source:
            self._data.append(item)

    def iterate(self, skip_empty=True):
        """ Iterate over data."""
        if skip_empty:
            for item in self._data:
                if not item:
                    continue
                yield item
        else:
            for item in self._data:
                yield item

    def iterate_with_func(self, pre_func, iter_func):
        """ Iterate over data with given iter_func.
        And data can be preprocessed with pre_func."""
        self._data = pre_func(self._data)
        for item in iter_func(self._data):
            yield item

    def do(self, cf, **args):
        """ Do something with data with given core function and args.
            And get a result of doing.
        """
        result = cf(self._data, **args)
        return result

    def process(self, cf, **args):
        """ Process data with given core function."""
        self._data = cf(self._data, **args)

    def clear(self):
        """ Remove data."""
        self._data = None

    def do_with_iter(self, cf, **args):
        """ Do something by iterating over data with given core function and args.
            And get a list of results of doing.
        """
        result = []
        for item in self._data:
            result.append(cf(item, **args))
        return result

    def process_with_iter(self, cf, **args):
        """ Process by iterating over data with given core function."""
        for i, item in enumerate(self._data):
            self._data[i] = cf(item, **args)

    def sort(self, sort_func, reverse=False):
        """ Sort data with sort_func and reversed param."""
        assert hasattr(sort_func, "__call__")
        self._data.sort(key=sort_func, reverse=reverse)

    @property
    def data(self):
        return self._data

    @property
    def N(self):
        return len(self._data)
   

def sc_sam_reader(sam_file):
    """
    """
    reader = SAMFileIO()
    