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
from trseeker.models.trf_model import TRModel
from trseeker.seqio.block_file import AbstractBlockFileIO
from trseeker.tools.sequence_tools import get_gc
from trseeker.seqio.mongodb_reader import MongoDBReader
from PyExp.readers.abstract_reader import sc_iter_filepath_folder
from trseeker.settings import load_settings

settings = load_settings()

class TRFFileIO(AbstractBlockFileIO):
    """  Working with raw ouput from TRF, where each block starts with '>' token.
    
    Public parameters:
    
    - self.use_mongodb -- Bool

    Public methods:
    
    - iter_parse(self, trf_file, filter=True)
    - parse_to_file(self, file_path, output_path, trf_id=0) -> trf_id
    
    Private methods:
    
    - _gen_data_line(self, data)
    - _filter_obj_set(self, obj_set)
    - _join_overlapped(self, obj1, obj2)
    
    Inherited public properties:
    
    - data  - iterable data, each item is tuple (head, body)
    - N     - a number of items in data
    
    Inherited public methods:
    
    - [OR] __init__(self)
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

    def __init__(self, use_mongodb=False):
        """ Overrided. Hardcoded start token."""
        token = "Sequence:"
        super(TRFFileIO, self).__init__(token)
        self.use_mongodb = use_mongodb
        if use_mongodb:
            db = MongoDBReader()
            db = db.get_trdb_conn()
            self.mongodb_table = db.tandem_repeats

    def iter_parse(self, trf_file, filter=True):
        """ Iterate over raw trf data and yield TRFObjs."""
        for head, body, start, next in self.read_online(trf_file):
            obj_set = []
            for line in self._gen_data_line(body):
                if not line:
                    continue
                trf_obj = TRModel()
                trf_obj.set_raw_trf(head, body, line)
                obj_set.append(trf_obj)
            if filter:
                # Unpack obj_set and sort by left, right induces
                obj_set = [(int(x.trf_l_ind), int(x.trf_r_ind), x) for x in obj_set]
                obj_set.sort()
                # Pack obj_set
                obj_set = [x[2] for x in obj_set]
                # Filter object set
                trf_obj_set = self._filter_obj_set(obj_set)
                obj_set = [x for x in trf_obj_set if x]
            yield obj_set

    def parse_to_file(self, file_path, output_path, trf_id=0, project=None):
        """ Parse trf file in tab delimited file."""
        if trf_id == 0:
            mode = "w"
        else:
            mode = "a"
        with open(output_path, mode) as fw:
            for trf_obj_set in self.iter_parse(file_path):
                for trf_obj in trf_obj_set:
                    trf_obj.trf_id = trf_id

                    if project:
                        trf_obj.set_project_data(project)

                    fw.write(str(trf_obj))

                    if self.use_mongodb:
                        self.write_to_mongodb(self.mongodb_table, trf_obj.__dict__)

                    trf_id += 1
        return trf_id


    def _gen_data_line(self, data):
        try:
            datas = re.compile('\n(\d.+?)\n', re.S).findall(data)
            for line in datas:
                yield line
        except:
            print "Failed parse data lines: %s" % (data)
            yield None

    def _filter_obj_set(self, obj_set):
        # Complex filter
        is_overlapping = False
        n = len(obj_set)

        obj_set.sort(key=lambda x: (x.trf_l_ind, x.trf_r_ind))
        for a in range(0, n):
            obj1 = obj_set[a]
            if not obj1:
                continue
            for b in range(a + 1, n):
                obj2 = obj_set[b]
                if not obj2:
                    continue
                # a ------ 
                # b ------
                if obj1.trf_l_ind == obj2.trf_l_ind and obj1.trf_r_ind == obj2.trf_r_ind:
                    # Check period
                    if obj1.trf_pmatch >= obj2.trf_pmatch:
                        obj_set[b] = None
                        continue
                    else:
                        obj_set[a] = None
                        continue
                # a ------ ------  ------- 
                # b ---       ---    ---
                if obj1.trf_l_ind <= obj2.trf_l_ind and obj1.trf_r_ind >= obj2.trf_r_ind:
                    obj_set[b] = None
                    continue
                # a ---       ---    ---
                # b ------ ------  ------- 
                if obj2.trf_l_ind <= obj1.trf_l_ind and obj2.trf_r_ind >= obj1.trf_r_ind:
                    obj_set[a] = None
                    continue
                # a ------ 
                # b    -----
                if obj1.trf_r_ind > obj2.trf_l_ind and obj1.trf_r_ind < obj2.trf_r_ind:
                    is_overlapping = True
                    continue
                # a ------ 
                # b         -----
                if obj1.trf_r_ind < obj2.trf_l_ind:
                    break
                # a         ------ 
                # b -----
                if obj2.trf_r_ind < obj1.trf_l_ind:
                    break
        obj_set = [a for a in obj_set if not a is None]
        n = len(obj_set)

        while is_overlapping:
            is_overlapping = False

            for a in range(0, n):
                obj1 = obj_set[a]
                if not obj1:
                    continue
                for b in range(a + 1, n):
                    obj2 = obj_set[b]
                    if not obj2:
                        continue
                    # a ------ 
                    # b         -----
                    if obj1.trf_r_ind < obj2.trf_l_ind:
                        break
                    # a         ------ 
                    # b -----
                    if obj2.trf_r_ind < obj1.trf_l_ind:
                        break
                    # a ------ 
                    # b    -----
                    if obj1.trf_r_ind > obj2.trf_l_ind and obj1.trf_r_ind < obj2.trf_r_ind:

                        overlap = float(abs(obj1.trf_r_ind - obj2.trf_l_ind))
                        min_length = min(obj1.trf_array_length, obj2.trf_array_length)
                        overlap_proc_diff = overlap * 1.0 / min_length
                        gc_dif = abs(obj1.trf_array_gc - obj2.trf_array_gc)

                        if overlap_proc_diff >= settings["trf_settings"]["overlapping_cutoff_proc"] \
                                    and gc_dif <= settings["trf_settings"]["overlapping_gc_diff"]:
                            is_overlapping = True
                            obj1 = self._join_overlapped(obj1, obj2)
                            obj2 = None
                        continue
                    # a ------ 
                    # b ------
                    if obj1.trf_l_ind == obj2.trf_l_ind and obj1.trf_r_ind == obj2.trf_r_ind:
                        # Check period
                        if obj1.trf_pmatch >= obj2.trf_pmatch:
                            obj_set[b] = None
                            continue
                        else:
                            obj_set[a] = None
                            continue
                    # a ------ ------  ------- 
                    # b ---       ---    ---
                    if obj1.trf_l_ind <= obj2.trf_l_ind and obj1.trf_r_ind >= obj2.trf_r_ind:
                        obj_set[b] = None
                        continue
                    # a ---       ---    ---
                    # b ------ ------  ------- 
                    if obj2.trf_l_ind <= obj1.trf_l_ind and obj2.trf_r_ind >= obj1.trf_r_ind:
                        obj_set[a] = None
                        continue

            obj_set = [a for a in obj_set if not a is None]
            n = len(obj_set)

        obj_set = [a for a in obj_set if not a is None]
        return obj_set

    def _join_overlapped(self, obj1, obj2):
        ''' Join overlapping sequences.'''
        obj1.trf_pmatch = int(obj1.trf_pmatch)
        obj2.trf_pmatch = int(obj2.trf_pmatch)
        obj1.trf_array_length = int(obj1.trf_array_length)
        obj2.trf_array_length = int(obj2.trf_array_length)

        obj1.trf_array = obj1.trf_array + obj2.trf_array[obj1.trf_r_ind - obj2.trf_l_ind:]

        if obj1.trf_array_length < obj2.trf_array_length:

            obj1.trf_consensus = obj2.trf_consensus
            obj1.trf_consensus_gc = obj2.trf_consensus_gc
            obj1.trf_period = obj2.trf_period
            obj1.trf_l_cons = obj2.trf_l_cons

            obj1.trf_n_a = obj2.trf_n_a
            obj1.trf_n_t = obj2.trf_n_t
            obj1.trf_n_c = obj2.trf_n_c
            obj1.trf_n_g = obj2.trf_n_g


        obj1.trf_pmatch = int((obj1.trf_pmatch * obj1.trf_array_length + obj2.trf_pmatch * obj2.trf_array_length) / (obj1.trf_array_length + obj2.trf_array_length))

        obj1.trf_l_ind = min(obj1.trf_l_ind, obj2.trf_l_ind)
        obj1.trf_r_ind = max(obj1.trf_r_ind, obj2.trf_r_ind)
        obj1.trf_array_length = obj1.trf_r_ind - obj1.trf_l_ind
        obj1.trf_n_copy = obj1.trf_array_length / obj1.trf_period

        obj1.trf_array_gc = get_gc(obj1.trf_array)

        obj1.trf_joined = 1

        return obj1

def sc_trf_mongodb_reader(project):
    """ Iter over trf file."""
    reader = TRFFileIO(use_mongodb=True)
    query = {"project":project["title"],
             "trf_array_length": {"$gt":project["parameters"]["array_length_for_large_trs"]}}
    for item in reader.read_from_mongodb(reader.mongodb_table, query):
        yield item


def sc_parse_raw_trf_folder(trf_raw_folder, output_trf_file, project=None):
    """ Parse raw TRF output in given folder to output_trf_file."""
    reader = TRFFileIO()
    trf_id = 1
    if os.path.isfile(output_trf_file):
        os.remove(output_trf_file)
    for file_path in sc_iter_filepath_folder(trf_raw_folder, mask="dat"):
        print "Start parse file %s..." % file_path
        trf_id = reader.parse_to_file(file_path, output_trf_file, trf_id=trf_id, project=project)

def sc_trf_mongodb_update(query_dict, update_dict):
    """ Shortcut for update TRF mongodb.
    
    Example:
    
    >>> query_dict = {'trf_array_length':96}
    >>> update_dict = {"id":0}
    
    """
    trf_file = TRFFileIO(use_mongodb=True)
    for x in trf_file.read_from_mongodb(trf_file.mongodb_table, query_dict):
        trf_file.update_mongodb(trf_file.mongodb_table, {"_id":x._id}, {"$set":update_dict})

