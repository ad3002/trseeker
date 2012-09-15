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

Organism model
--------------

::

	from trseeker.models.organism_model import OrganismModel		

Attributes:

- organism_taxon
- organism_common_name
- organism_acronym
- organism_description
- organism_wgs_projects
- organism_genome_assemblies

Dataset model
-------------

::

	from trseeker.models.dataset_model import DatasetModel		

Attributes:

- dataset_taxon
- dataset_id
- dataset_sources
- dataset_description
- dataset_gc (float)
- dataset_length (int)
- dataset_trs_n (int)
- dataset_trs_length (int)
- dataset_trs_mean_gc (float)
- dataset_trs_fraq (float)

Blast Results Model
-------------------

::

	from trseeker.models.blast_model import BlastResultModel		

Attributes:

- query_id (int)
- query_gi (int)
- query_ref
- subject_id
- subject_gi(int)
- subject_ref
- query_start (int)
- query_end (int)
- subject_start (int)
- subject_end (int)
- evalue  (float)
- bit_score (flaot)
- score (int)
- alignment_length (int)
- proc_identity (float)
- identical (int)
- mismatches (int)
- positives (int)
- gap_opens (int)
- gaps (int)
- proc_positives (float)
- frames
- query_frame (int)
- subject_frame (int)
- fraction_of_query (float)  

Additional functions:

- read_blast_file(blast_file, length), return subject_ref -> list of matches (BlastResultModel models).

::

	from trseeker.models.blast_model import read_blast_file

	ref_to_blast_obj = read_blast_file(file_name)

Chromosome model
----------------

::

	from trseeker.models.chromosome_model import ChromosomeModel

Attributes:

- chr_genome
- chr_number
- chr_taxon
- chr_prefix
- chr_gpid
- chr_acronym
- chr_contigs
- chr_length
- chr_mean_gc
- chr_trs_all
- chr_trs_3000
- chr_trs_all_proc
- chr_trs_3000_proc
- chr_trs_all_length
- chr_trs_3000_length
- genome_gaps
- chr_sum_gc
	
WGS model
---------

::

	from trseeker.models.wgs_model import WGSModel

Attributes:

- wgs_prefix
- wgs_taxon
- wgs_gpid
- wgs_acronym
- wgs_contigs (int)
- wgs_length (int)
- wgs_mean_gc (float)
- wgs_trs_all (int)
- wgs_trs_3000 (int)
- wgs_trs_1000 (int)
- wgs_trs_500 (int)
- wgs_trs_all_proc (float)
- wgs_trs_3000_proc (float)
- wgs_trs_1000_proc (float)
- wgs_trs_500_proc (float)
- wgs_trs_all_length (int)
- wgs_trs_3000_length (int)
- wgs_trs_1000_length (int)
- wgs_trs_500_length (int)
- wgs_sum_gc (float)

Methods:

- wgs-obj.clear_trf(), clear trf information (set to 0)

Genome model
------------

::

	from trseeker.models.genome_model import GenomeModel

- genome_taxon
- genome_prefix
- genome_gpid
- genome_acronym
- genome_chromosomes
- genome_contigs
- genome_length
- genome_mean_gc
- genome_trs_all
- genome_trs_3000
- genome_trs_all_proc
- genome_trs_3000_proc
- genome_trs_all_length
- genome_trs_3000_length
- genome_gaps
- genome_sum_gc

Ngram model
-----------

::

	from trseeker.models.ngram_model import NgramModel

	ngram_obj = NgramModel(seq_f, seq_r)
	ngram_obj.add_tr(trf_obj, tf)

	print ngram_obj
	>>> '[fseq]\t[rseq]\t[tf]\t[df]\t[len taxons]\t[len fams]\n'

	print ngram_obj.get_families()
	>>> ???

Attributes

- seq_r
- seq_f
- tf (float)
- df (int)
- taxons (set)
- trs (set)
- families (dict)

Ngrams model
------------

::

	from trseeker.models.ngrams_model import NgramModel

Attributes

- id (int)
- rev_id (int)
- ngram
- rev_ngram
- tf (float)
- df (int)
- etf (float)
- edf (int)

NgramToTRModel model
--------------------

::

	from trseeker.models.ngrams_model import NgramToTRModel

Attributes

- id (int)
- trids (list int)
- tfs (list int)

Additional function

- sc_ngram_reader(file_name), yield NgramModel
- sc_ngram_trid_reader(file_name), yield (ngram id, [(seq id, tf), ...])


IO functions
============

Toolkit
=======
