#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 05.06.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com
"""
Classes:
    
- NCBIFtpIO(AbstractFtpIO)
    
"""
from trseeker.seqio.ftp_io import AbstractFtpIO
import os

class NCBIFtpIO(AbstractFtpIO):
    """ Class for working with data via FTP.
        
    Public methods:
    
    - download_wgs_fasta(self, wgs_list, file_suffix, output_folder)
    
    Ingerited public methods:
    
    - [OR] __init__(self)
    - connect(self)
    - cd(self, list_of_dirs)
    - ls(self)
    - get(self, file, output_file)
    - unzip(self, file_name)
    - download_all_wgs_in_fasta(self, output_folder)
    - download_all_wgs_in_gbff(self, output_folder)
    """

    def __init__(self):
        """ Set ftp parameters.
        
        Keyword arguments:
        
        - ftp_address   -- address of ftp server
        """
        ftp_address = "ftp.ncbi.nlm.nih.gov"
        super(NCBIFtpIO, self).__init__(ftp_address)

    def download_wgs_fasta(self, wgs_list, file_suffix, output_folder):
        """ Download and unzip (gz) wgs file from wgsfor given wgs_list projects and file suffix.
            
        Parameters:

        - wgs_list  -- list of WGS projects, eg ["AABB",]
        - file_suffix - extension for file_format
        - output_folder -- folder for fasta files
         """
        
        print "WGS list: ", wgs_list
        print "File suffix: ", file_suffix

        print "Read NCBI's ftp..."
        path = ["genbank", "wgs"]
        self.connect()
        self.cd(path)
        files = self.ls()
        files = [ item for item in files if [True for project in wgs_list if project in item] ]
        files = [ item for item in files if file_suffix in item ]

        print "Files to download: ", files

        for file_name in files:
            output_file = os.path.join(output_folder, file_name)
            print "Start download: %s ..." % file_name,
            print " to %s" % output_file
            self.get(file_name, output_file)
            print "Unzip..."
            self.unzip(output_file)
            os.unlink(output_file)

    def download_all_wgs_in_fasta(self, output_folder):
        """ Download all WGS files from NCBI in fasta format."""
        file_suffix = "fsa_nt"
        path = ["genbank", "wgs"]
        self.connect()
        self.cd(path)
        files = self.ls()
        files = [ item for item in files if file_suffix in item ]
        n = len(files)
        for i, file in enumerate(files):
            output_file = os.path.join(output_folder, file)
            print "Start download: %s ..." % file,
            print " to %s" % output_file

            if os.path.isfile(output_file):
                print "--> was downloaded early"
                continue
            self.get(file, output_file)
            print "Unzip..."
            self.unzip(output_file)
            print "Done %i from %s" % (i, n)

    def download_all_wgs_in_gbff(self, output_folder):
        """ Download all WGS files from NCBI in genbank format."""
        file_suffix = "gbff"
        path = ["genbank", "wgs"]
        self.connect()
        self.cd(path)
        files = self.ls()
        files = [ item for item in files if file_suffix in item ]
        for file in files:
            output_file = os.path.join(output_folder, file)
            print "Start download: %s ..." % file,
            print " to %s" % output_file
            self.get(file, output_file)
            print "Unzip..."
            self.unzip(output_file)

    def download_trace_data(self, taxon, output_folder):
        """ Download Trace data."""
        file_suffix = "gbff"
        path = ["pub", "TraceDB", taxon]
        self.connect()
        self.cd(path)
        files = self.ls()
        files = [ item for item in files if "fasta" in item or "clip" in item ] 
        index = set()
        for file in files:
            output_file = os.path.join(output_folder, file)
            print "Start download: %s ..." % file,
            print " to %s" % output_file
            self.get(file, output_file)
            print "Unzip..."
            self.unzip(output_file)
            index.add(file.split(".")[-2])
        return index
