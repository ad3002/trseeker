#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 05.06.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
SRA files.
TODO: implement this.
'''

class FastqObj(object):

	def __init__(self, head, seq, strain, qual):
		self.head = head
		self.seq = seq.lower()
		self.qual = qual
		self.strain = strain
		self.id = head.strip().split()[0].replace("@","")

	@property
	def fastq(self):
		return "%s%s%s%s" % (
					self.head,
					self.seq,
					self.strain,
					self.qual,
			)

	@property
	def sequence(self):
		return self.seq.strip().lower()

	@property
	def fasta(self):
		return ">%s\n%s\n" % (
					self.head.strip(),
					self.seq.strip(),
					)

def fastq_reader(fastq_file):
	with open(fastq_file) as fh:
		while True:
			try:
				head = fh.readline()
				seq = fh.readline()
				strain = fh.readline()
				qual = fh.readline()
				fastq_obj = FastqObj(head, seq, strain, qual)
				yield fastq_obj
			except:
				break