===========================================================
Trseeker: a Python framework for working with noncoding DNA
===========================================================

Framework состоит из следующих частей

- модели данных
- чтение/запись данных согласно моделям
- инструменты работы с этими данными

Настройки фреймворка
--------------------

В файле settings.py

::

	SETTINGS_FILENAME = "settings.yaml"
	NGRAM_LENGTH = 12
	NGRAM_N = 100000000

Настройки можно прочитать и записать:

::
	
	from trseeker.settings import load_settings
	from trseeker.settings import save_settings

	settings_dict = load_settings()
	save_settings(settings_dict)

Настройки содержат следующие параметры:
	
::

	trseeker:
	    os: linux64
	    root_dir: /root/Dropbox/workspace/trseeker
	    work_dir: /home
	blast_settings:
	    blast_location: 
	    repbase_db_folder: /home/rebase_blast_db/repbase
	    blast_e: 0.000000000001
	    blast_b: 1000000
	    blast_v: 1000000
	trf_settings:
	    trf_location: /root/trf404.linux64
	    trf_match: 2
	    trf_mismatch: 5
	    trf_indel: 7
	    trf_p: 80
	    trf_q: 10
	    trf_threshold: 50
	    trf_length: 2000
	    trf_param_postfix: 2.5.7.80.10.50.2000
	    trf_masked_file: False
	    trf_flanked_data: True
	    trf_data_file: True
	    trf_nohtml: True
	    overlapping_cutoff: 10
	    overlapping_cutoff_proc: 30
	    overlapping_gc_diff: 0.05
	ngrams_settings:
	    ngram_length: 12
	    ngram_m: 10000000

Avaliable Models
================

DNA Sequence
------------

::

	from trseeker.models.sequence_model import SequenceModel

Attribites

- seq_gi (int)
- seq_ref
- seq_description
- seq_sequence
- seq_length (int)
- seq_gc (float)
- seq_revcom, reverce complement
- seq_gapped (int)
- seq_chr
- seq_head
- seq_start_position (int)
- seq_end_position (int)

Properties

- length (self.seq_length)
- sequence (self.seq_sequence)
- fasta 

::

	print seq_obj.fasta
	>>> ">[seq_ref]\n[seq_sequence]\n"

- sa_input

::

	print seq_obj.sa_input
	>>> "[seq_sequence]$"

- ncbi_fasta

::

	print seq_obj.ncbi_fasta
	>>> ">gi|[seq_gi]|ref|[seq_ref]|[seq_description]\n[seq_sequence]\n"

Methods

- add_sequence_revcom()
- set_dna_sequence(self, title, sequence, description=None)

::

	self.seq_ref = title
    
- set_ncbi_sequence(self, head, sequence)

::
	
	(self.seq_gi, self.seq_ref, self.seq_description) = parse_fasta_head(head)

Chromosome name is **?** or setted with parse_chromosome_name(head).

- set_gbff_sequence(self, head, sequence)

Head is a dictionary with gi, ref, description keys.

Chromosome name is **?** or setted with parse_chromosome_name(head["description"]).

Sequence is cleared with clear_sequence(s) function. Lowercase and all non-DNA characters replacing with **n**. If the sequence has **n** then it is gapped.

TRF results
-----------

::

	from trseeker.models.trf_model import TRModel

Attributes:

- project, project name
- id (int)
- trf_id (int)
- trf_l_ind (int)
- trf_r_ind (int)
- trf_period (int)
- trf_n_copy (float)
- trf_pmatch (float)
- trf_pvar (float)
- trf_consensus
- trf_array
- trf_array_gc (float)
- trf_consensus_gc (float)
- trf_gi
- trf_head
- trf_param
- trf_array_length (int)
- trf_chr
- trf_joined (int)
- trf_superfamily
- trf_superfamily_self
- trF_superfamily_ref
- trf_family
- trf_subfamily
- trf_subsubfamily
- trf_family_network
- trf_family_self
- trf_family_ref
- trf_hor (int)
- trf_n_chrun (int)
- trf_chr_refgenome
- trf_bands_refgenome
- trf_repbase
- trf_strand

Methods

- set_project_data(project), set self.project to given project
- set_raw_trf(head, body, line), head, body and line from TRF parser
- get_index_repr()

::
	
	print trf_obj.get_index_repr()
	'''
	Tab delimted string with \n-symbol:
	trf_id
	trf_period
	trf_array_length
	trf_pvar
	trf_gi
	trf_l_ind
	trf_r_ind
	trf_chr
	'''

- get_numerical_repr()

::
	print trf_obj.get_numerical_repr()
	>>> [trf_period]\t[trf_array_length]\t[trf_array_gc]\n

- get_fasta_repr(), where head is trf_obj.trf_id and sequence is trf_obj.trf_array
- get_monomer_fasta_repr(), where head is trf_obj.trf_id and sequence is trf_obj.trf_consensus
- get_family_repr()

::
	
	print trf_obj.get_family_repr()
	'''
	Tab delimted string with \n-symbol:
	trf_id
	trf_period
	trf_array_length
	trf_array_gc
	trf_pvar
	trf_gi
	trf_l_ind
	trf_r_ind
	trf_chr
	trf_repbase
	trf_superfamily
	trf_family
	trf_subfamily
	'''

For network slice added one more index - gid (group id)

::

	from trseeker.models.trf_model import NetworkSliceModel	

	slice_obj = NetworkSliceModel()



IO functions
============

Toolkit
=======
