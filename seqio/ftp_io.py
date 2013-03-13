#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 05.06.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
"""
Classes:
    
- AbstractFtpIO(object)
    
"""
import ftplib
import os
import gzip

class AbstractFtpIO(object):
    """ Class for working with data via FTP.
        
    Public methods:
    
    - __init__(self, ftp_address=None)
    - connect(self)
    - cd(self, list_of_dirs)
    - ls(self)
    - get(self, file, output_file)
    - unzip(self, file_name)
        
    """

    def __init__(self, ftp_address=None):
        """ Set ftp parameters.
        
        Keyword arguments:
        
        - ftp_address   -- address of ftp server
        """
        self.ftp_address = ftp_address
        self.ftp = None

    def connect(self):
        """ Connect anonymously to server."""
        print "Connect to", self.ftp_address
        self.ftp = ftplib.FTP(self.ftp_address)
        self.ftp.login()

    def cd(self, list_of_dirs):
        """ Go to dir by given path (list of dirs)."""
        for dirname in list_of_dirs:
            self.ftp.cwd(dirname)

    def ls(self):
        """ Return list of files filtered with filter_func in folder.
        """
        return self.ftp.nlst()

    def get(self, file, output_file):
        """ Save file in output file."""
        with open(output_file, "ab") as fh:
            self.ftp.retrbinary("RETR %s" % file, lambda x: fh.write(x))

    def unzip(self, file_name):
        ''' Unzip file.'''
        if not file_name.endswith(".gz"):
            return
        fh = gzip.GzipFile(fileobj=open(file_name, 'rb'))
        fw = file(file_name.replace('.gz', ''), 'wb')
        for line in fh:
            fw.write(line)
        fh.close()
        fw.close()


