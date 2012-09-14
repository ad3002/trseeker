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
~~~~~~~~~~

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
~~~~~~~~~~

- length (self.seq_length)
- sequence (self.seq_sequence)
- fasta 

::

	print seq_obj.fasta
	>>> ">[seq_ref]\n[seq_sequence]\n"




IO functions
============

Toolkit
=======
