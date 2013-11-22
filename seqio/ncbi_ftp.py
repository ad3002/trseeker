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
from trseeker.seqio.ftp_io import AbstractFtpIO, download_with_aspera_from_ncbi
import os
import re

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

    def __init__(self, aspera=True):
        """ Set ftp parameters.
        
        Keyword arguments:
        
        - ftp_address   -- address of ftp server
        """
        ftp_address = "ftp.ncbi.nlm.nih.gov"
        super(NCBIFtpIO, self).__init__(ftp_address)
        self.aspera = aspera

    def download_wgs_fasta(self, wgs_list, file_suffix, output_folder, unzip=False):
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
            if self.aspera:
                ftp_address = self.ftp_address + "/genbank/wgs/" + file_name
                self.download_with_aspera(output_file, "NCBI", ftp_address)
            else:
                print "Start download: %s ..." % file_name,
                print " to %s" % output_file
                self.get(file_name, output_file)
                if unzip:
                    print "Unzip..."
                    self.unzip(output_file)
                    os.unlink(output_file)

    def download_all_wgs_in_fasta(self, output_folder, aspera=False, get_size=False, unzip=False):
        """ Download all WGS files from NCBI in fasta format."""
        file_suffix = "fsa_nt"
        path = ["genbank", "wgs"]
        self.connect()
        self.cd(path)
        files = self.ls()
        files = [ item for item in files if file_suffix in item ]
        n = len(files)
        total_size = 0
        for i, file in enumerate(files):
            if get_size:
                total_size += self.size(file)
            print i, total_size / 8589934592, "\r"
            output_file = os.path.join(output_folder, file)
            if self.aspera:
                ftp_address = self.ftp_address + "/genbank/wgs/" + file
                self.download_with_aspera(output_file, "NCBI", ftp_address)
            else:
                print "Start download: %s ..." % file,
                print " to %s" % output_file
                if os.path.isfile(output_file):
                    print "--> was downloaded early"
                    continue
                if aspera:
                    download_with_aspera_from_ncbi("/genbank/wgs/"+file, output_file)
                else:
                    self.get(file, output_file)
                if unzip:
                    print "Unzip..."
                    self.unzip(output_file)
                    print "Done %i from %s" % (i, n)

    def download_chromosomes_in_fasta(self, ftp_address, name, output_folder, unzip=False):
        """ Download all WGS files from NCBI in fasta format."""

        if self.aspera:
            self.download_with_aspera(output_folder, "NCBI", ftp_address)
            return

        paths = ftp_address.split("/")
        if paths[-1].endswith("/"):
            self.ftp_address, folders = paths[0], paths[1:]
        else:
            self.ftp_address, folders = paths[0], paths[1:-1]
        self.connect()
        self.cd(folders)
        files = self.ls()
        if name:
            files = [ item for item in files if re.search(name, item)]
        n = len(files)
        for i, file in enumerate(files):
            output_file = os.path.join(output_folder, file)
            print "Start download: %s ..." % file,
            print " to %s" % output_file

            if os.path.isfile(output_file):
                print "--> was downloaded early"
                continue
            self.get(file, output_file)
            if unzip:
                print "Unzip..."
                self.unzip(output_file)
            print "Done %i from %s" % (i, n)

    def download_all_wgs_in_gbff(self, output_folder, unzip=False):
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
            if unzip:
                print "Unzip..."
                self.unzip(output_file)

    def download_trace_data(self, trace_folder, output_folder, unzip=False):
        """ Download Trace data."""
        file_suffix = "gbff"
        path = ["pub", "TraceDB", trace_folder]
        self.connect()
        self.cd(path)
        files = self.ls()
        files = [ item for item in files if "fasta" in item or "clip" in item ] 
        index = set()
        for file_name in files:
            if self.aspera:
                ftp_address = self.ftp_address + "/pub/TraceDB/%s/%s" % (trace_folder, file_name)
                self.download_with_aspera(output_folder, "NCBI", ftp_address)
            else:
                print "Start download: %s ..." % file_name,
                print " to %s" % output_file
                self.get(file_name, output_file)
                if unzip:
                    print "Unzip..."
                    self.unzip(output_file)
                index.add(file_name.split(".")[-2])
        return index

    def download_sra_from_ddbj(self, ftp_address, output_folder):
        '''
        '''
        file_name = ftp_address.rsplit("/")[-1]
        if self.aspera:
            self.download_with_aspera(output_folder, "NCBI", ftp_address)
        else:
            output_file = os.path.join(output_folder, file_name)
            paths = ftp_address.split("/")
            self.ftp_address, folders = paths[0], paths[1:]
            self.connect()
            self.cd(folders)
            print "Start download as is: %s ..." % file_name,
            print " to %s" % output_file
            self.get(file_name, output_file)

    def download_with_aspera(self, local_path, remote_server, remote_path):
        '''
        TODO: move path to settings
        '''
        if remote_server == "NCBI":
            remote_server = "anonftp@ftp-private.ncbi.nlm.nih.gov"
            remote_path = "/%s" % "/".join(remote_path.split("/")[1:])
        params = {
            "location": "/home/akomissarov/.aspera/connect/bin/ascp",
            "putty_keys": "/home/akomissarov/.aspera/connect/etc/asperaweb_id_dsa.putty",
            "remote_server": remote_server,
            "remote_path": remote_path,
            "local_path": local_path,
        }
        command = "%(location)s -QT -l640M -d --overwrite=diff -i %(putty_keys)s %(remote_server)s:%(remote_path)s %(local_path)s" % params
        print command
        os.system(command)