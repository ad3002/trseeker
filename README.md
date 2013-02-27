# Trseeker: a Python framework for working with noncoding DNA

Framework состоит из следующих частей:

- модели данных
- чтение/запись данных согласно моделям
- инструменты работы с этими данными

## Настройки фреймворка

В файле settings.py:

```python
SETTINGS_FILENAME = "settings.yaml"
NGRAM_LENGTH = 23
NGRAM_N = 100000000
```

Настройки можно прочитать и записать:

```python
from trseeker.settings import load_settings
from trseeker.settings import save_settings

settings_dict = load_settings()
save_settings(settings_dict)
```

Настройки содержат следующие параметры:

```yaml	
trseeker:
    os: linux64
    root_dir: /root/Dropbox/workspace/trseeker
    work_dir: /home
blast_settings:
    blast_location: 
    repbase_db_folder: /home/rebase_blast_db/repbase
    blast_e: 1e-20
    blast_b: 20000000
    blast_v: 20000000
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
    ngram_length: 23
    ngram_m: 10000000
```

## Avaliable Models

### DNA Sequence

```python
from trseeker.models.sequence_model import SequenceModel
```

Attribites:

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

```python
print seq_obj.fasta
>>> ">[seq_ref]\n[seq_sequence]\n"
```

- sa_input

```python
print seq_obj.sa_input
>>> "[seq_sequence]$"
```

- ncbi_fasta

```python
print seq_obj.ncbi_fasta
>>> ">gi|[seq_gi]|ref|[seq_ref]|[seq_description]\n[seq_sequence]\n"
```

Methods:

- add_sequence_revcom()
- set_dna_sequence(self, title, sequence, description=None)

```python
self.seq_ref = title
```

- set_ncbi_sequence(self, head, sequence)

Chromosome name is **?** or setted with parse_chromosome_name(head).

```python	
(self.seq_gi, self.seq_ref, self.seq_description) = parse_fasta_head(head)
```

- set_gbff_sequence(self, head, sequence)

Head is a dictionary with gi, ref, description keys.

Chromosome name is **?** or setted with parse_chromosome_name(head["description"]).

Sequence is cleared with clear_sequence(s) function. Lowercase and all non-DNA characters replacing with **n**. If the sequence has **n** then it is gapped.

### TRF results

```python
from trseeker.models.trf_model import TRModel
```

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
- trf_entropy (float)
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

```python	
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
```

- get_numerical_repr()

```python
print trf_obj.get_numerical_repr()
>>> [trf_period]\t[trf_array_length]\t[trf_array_gc]\n
```

- get_fasta_repr(), where head is trf_obj.trf_id and sequence is trf_obj.trf_array
- get_monomer_fasta_repr(), where head is trf_obj.trf_id and sequence is trf_obj.trf_consensus
- get_family_repr()

```python
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
```

For network slice added one more index - gid (group id)

```python
from trseeker.models.trf_model import NetworkSliceModel	

slice_obj = NetworkSliceModel()
```

### Organism model

```python
from trseeker.models.organism_model import OrganismModel		
```

Attributes:

- organism_taxon
- organism_common_name
- organism_acronym
- organism_description
- organism_wgs_projects
- organism_genome_assemblies

### Dataset model

```python
from trseeker.models.dataset_model import DatasetModel		
```

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

### Blast Results Model

```python
from trseeker.models.blast_model import BlastResultModel		
```

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

```python
from trseeker.models.blast_model import read_blast_file

ref_to_blast_obj = read_blast_file(file_name)
```

### Chromosome model

```python
from trseeker.models.chromosome_model import ChromosomeModel
```

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
	
### WGS model

```python
from trseeker.models.wgs_model import WGSModel
```

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

Clear trf information (set to 0):

```python
wgs_obj.clear_trf()
```

### Genome model

```python
from trseeker.models.genome_model import GenomeModel
```

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

### Ngram model

```python
	from trseeker.models.ngram_model import NgramModel
```

```python
ngram_obj = NgramModel(seq_f, seq_r)
ngram_obj.add_tr(trf_obj, tf)

print ngram_obj
>>> '[fseq]\t[rseq]\t[tf]\t[df]\t[len taxons]\t[len fams]\n'

print ngram_obj.get_families()
>>> ???
```

Attributes

- seq_r
- seq_f
- tf (float)
- df (int)
- taxons (set)
- trs (set)
- families (dict)

### Ngrams model

```python
from trseeker.models.ngrams_model import NgramModel
```

Attributes

- id (int)
- rev_id (int)
- ngram
- rev_ngram
- tf (float)
- df (int)
- etf (float)
- edf (int)

### NgramToTRModel model

```python
from trseeker.models.ngrams_model import NgramToTRModel
```

#### Attributes

- id (int)
- trids (list int)
- tfs (list int)

#### Additional function

Yield NgramModel:

```python
sc_ngram_reader(file_name)
```

Yield (ngram id, [(seq id, tf), ...]):

```python
sc_ngram_trid_reader(file_name)
```

## IO functions

### Tab file

```python
from trseeker.seqio.tab_file import TabDelimitedFileIO

reader = TabDelimitedFileIO(skip_first=False, format_func=None, delimiter='\t', skip_startswith='#')

reader.read_from_file(file_name)
```

Avaliable all functions from parent class AbstractFileIO from PyExp package.

Useful functions:

- sc_iter_tab_file(input_file, data_type, remove_starts_with=None)
- sc_iter_simple_tab_file(input_file)
- sc_read_dictionary(dict_file, value_func=None)
- sc_read_simple_tab_file(input_file)

```python	
from trseeker.seqio.tab_file import sc_iter_tab_file

for wgs_obj in sc_iter_tab_file(file_name, WGSModel):
	print wgs_obj.wgs_taxon
```

```python
from trseeker.seqio.tab_file import sc_iter_simple_tab_file

for item in sc_iter_simple_tab_file(file_name):
	print item[0]
```

```python
from trseeker.seqio.tab_file import sc_read_dictionary

for data in sc_read_dictionary(file_name):
	for key, value in data.items():
		print key, value
```

```python
from trseeker.seqio.tab_file import sc_read_simple_tab_file

for data in sc_read_simple_tab_file(file_name):
	for item in data:
		print "id = ", item[0]	
```

### Block file

```python
from trseeker.seqio.block_file import AbstractBlockFileIO

reader = AbstractBlockFileIO(token, **args)
for (head, body, start, next) in reader.read_online(file_name):
	pirnt head
```

Avaliable all functions from parent class AbstractFileIO from PyExp package.

### Fasta file

```python
from trseeker.seqio.fasta_file import FastaFileIO

reader = FastaFileIO()
```

#### Useful functions:

- sc_iter_fasta(file_name)
- sc_iter_fasta_simple(file_name)
- save_all_seq_with_exact_substring(fasta_file, substring, output_file)

```python
from trseeker.seqio.fasta_file import sc_iter_fasta

for seq_obj in sc_iter_fasta(file_name):
	print seq_obj.sequence
```

```python
from trseeker.seqio.fasta_file import sc_iter_fasta_simple

for (gi, sequence) in sc_iter_fasta_simple(file_name):
	print gi
```

```python
from trseeker.seqio.fasta_file import save_all_seq_with_exact_substring

save_all_seq_with_exact_substring(fasta_file, substring, output_file)
```

### FTP IO

```python	
from trseeker.seqio.frt_io import AbstractFtpIO

reader = AbstractFtpIO(ftp_address=address)

reader.connect()

reader.cd(['/home', 'user', 'data'])

reader.ls()
>>> ['readme.txt', 'data.fa']

reader.get(file, output_file)

reader.unzip(file_name)
```

### NCBI ftp

```python
from trseeker.seqio.ncbi_ftp import NCBIFtpIO

reader = NCBIFtpIO()

reader.download_wgs_fasta(wgs_list, file_suffix, output_folder)

reader.download_all_wgs_in_fasta(output_folder)

reader.download_all_wgs_in_gbff(output_folder)

reader.download_trace_data(taxon, output_folder)
```

### Mongo db reader

Not implemented yet.

```python	
from trseeker.seqio.mongodb_reader import MongoDBReader

reader = MongoDBReader()
db = reader.get_trdb_conn()
```

### GBFF file

```python
from trseeker.seqio.gbff_file import GbffFileIO

reader = GbffFileIO()
```

#### Useful functions:

- sc_iter_gbff(file_name)
- sc_iter_gbff_simple(file_name)
- sc_parse_gbff_in_folder(gbbf_folder, fasta_folder, fasta_postfix, mask)

```python
from trseeker.seqio.gbff_file import sc_iter_gbff

for seq_obj in sc_iter_gbff(file_name):
	print seq_obj.length
```

```python
from trseeker.seqio.gbff_file import sc_iter_gbff_simple

for (gi, sequence) in sc_iter_gbff_simple(file_name):
	print gi
```

```python
from trseeker.seqio.gbff_file import sc_parse_gbff_in_folder

sc_parse_gbff_in_folder("/home/user/", "/home/user/fasta", "fa", "mouse")
```

### TRF file

#### Useful functions:

- sc_parse_raw_trf_folder(trf_raw_folder, output_trf_file, project=None)

```python
from trseeker.seqio.trf_file import TRFFileIO

from trseeker.seqio.trf_file import sc_parse_raw_trf_folder

for trf_obj in TRFFileIO(file_name, filter=True):
	print trf_obj.trf_id
```

```python
sc_parse_raw_trf_folder(trf_raw_folder, output_trf_file, project="mouse_genome")
```

При чтение данных TRF происходит их фильтрация по следующим параметрам:

1. Убираются все вложенные поля меньшей длины.
2. Если поля overlapping, то сейчас ничего не делается, а раньше если 

```python
overlap_proc_diff >= settings["trf_settings"]["overlapping_cutoff_proc"] 
gc_dif <= settings["trf_settings"]["overlapping_gc_diff"]
``` 
поля объяединяются в одно поле. Иначе считаем, что это отдельные поля.
3. Если поля совпадают, то выбирается то которое с большим trf_pmatch.

TODO: убрать пересечение, так как это видимо в том числе и ошибки ассмеблера и это будет мешать корректной классификации. Убрать в отдельную часть.

### TR file

```python
from trseeker.seqio.tr_file import *
```

Functions:

- save_trs_dataset(trs_dataset, output_file)
- get_trfid_obj_dict(trf_large_file)

```python
trid2trfobj = get_trfid_obj_dict(trf_large_file)
```

- get_all_trf_objs(trf_large_file)

```python
trf_obj_list = get_all_trf_objs(trf_large_file)
```

- read_trid2meta(file_name)

Load trid to full index dictionary as string:

```python
read_trid2ngrams(annotation_ngram_folder, trf_large_file)
```

TODO: rewrite this with AbstractReaders

### Ngram file

TODO: rewrite this with AbstractReaders

```python
from trseeker.seqio.ngram_file import save_ngram_index

save_ngram_index(ngram_index_file,
		       hash2id,
		       result_tf,
                     	       result_df, 
                     	       result_rtf, 
                     	       result_rdf,
	                      seen_rev, 
	                      hash2rev)

save_ngram_pos_index(ngram_trids_file, id2trids, id2trid2tf)    

save_distance_data(dist_file, distances)
```

### Blast results file

```python	
from trseeker.seqio.blast_file import get_blast_result

get_blast_result(blast_file, length, gap_size=1000, min_align=500, min_length=2400, format_function=None)
```

Dataset updating functions:

```python
update_with_repbase_blast_result(trs_dataset, annotation_self_folder, filters)

# inside function
# result is semicolon-delimited
trs_dataset[i].trf_repbase = result

# filters example
filters = {
	"blast_gap_size": 300,
	"min_align": 1000,
	"min_length": 3000,
}
```

```python
update_with_self_blast_result(trs_dataset, annotation_self_folder, _get_filters)
	
# inside function
# result comma-delimited
# _get_filters(array_length) generate parameters by array_length
trs_dataset[i].trf_family_self = result
```

```python
update_with_ref_blast_result(trs_dataset, annotation_self_folder, filters)

# inside function
# comma-delimited data
trs_dataset[i].trf_family_ref = result
```

TODO: refractor to more generic form with setattr(...)

Self and ref annotation comma-delimited. Repbase annotation semicolon-delimited. Avalibale format functions:

- format_function_repbase(blast_dataset)
- format_function_self(blast_dataset) [DEFAULT]

While blast results parsing:
- skipped all inner matches
- joined overlapped
- skipped all alignments with length less thatn min_align paramter
- joined gapped if gap less than gap_size paramter

### Sra file

```python
from trseeker.seqio.sra_file import FastqObj

fastq_obj = FastqObj(head, seq, srain, qual_str)
print fastq_obj.fastq
```

Additional functions:

- fastq_reader

```python
for fastq_obj in fastq_reader(fastq_file):
	print fastq_obj.seq
```

## Toolkit

### Tulip files

Format and write network data to tulip format from distance_file with given cutoff (defaul=90).
    
- distance_file: (node_a, node_b, distance)
- tulip_file_output: network in tulip format
- cutoff: distance cutoff

```python
from trseeker.tools.tulip_tools import *

format_distances_to_tulip(distance_file, tulip_file_output, cutoff=90)
```

### TRs groups

```python
from trseeker.tools.trs_groups import *
```python

Get next family index. Index limited to latin alphabet characters. Otherwise will return X1,X2 and so on

```python
letter = get_index(i)
```python

Get most frequent element of list.

```python
element = get_popular(s)
```

Get unit and letter for family. Unit picked by most common pmatch value. A seen_units is a dictionary unit_size to frequence, it is needed for letter choosing.

```python
unit, letter, seen_units = get_family_name(ids, seen_units)	
```

Join families with common members.

```python
join_families_with_common(families)
```

### Sequence tools

```python
from trseeker.tools.sequence_tools import *
```

Return complementary sequence.

```python
revcom = get_revcomp(sequence)
```

Return normalized sequence with following rules:
1) if T > A then return reverse complement
2) if A == T and G > C  then return reverse complement
3) else return sequence

```python
sequence = fix_strand(sequence)
```

Count GC content.
```python
gc = get_gc(sequence)
```

Check for N|n in sequence. Return n(with N) or w(whole).
```python
bool_value = check_gapped(sequence)
```

Return subsequence.
```python
get_subseq(seq, start, end)
```

Clear sequence (ACTGN alphabet)
1) lower case
2) remove any gaps
3) remove all letters except atgcn

```python
sequence = clear_sequence(sequence)
```

Return list of sequences splited by pattern with added end. Pattern is regexp.
```python
restriction(sequence, pattern, end="")
```
Return number of fragments after restriction of sequence with given pattern. Pattern is regexp.
```python
restriction_fragments_n(sequence, pattern)
```
Check tandem repeat synonims  between two sequence with same length.
```python
check_cyclic_repeats(seq_a, seq_b)
```
Return sequence with n mutations.
```python
random_mutation(seq, n, alphabet="actgn +")
```
Return consensus string for given list of strings.
```python
get_consensus(strs)
```

### Various useful functions

	from trseeker.tools.seqfile import *

- save_list(file_name, data)
- save_dict(file_name, dictionary)
- save_sorted_dict(file, d, by_value=True, reverse=True)
- count_lines(file)
- sort_file_by_int_field(file_name, field)

### Other useful functions (other_tools.py)

	from trseeker.tools.other_tools import *

- sort_dictionary_by_key(d, reverse=False)
Sort dictionary by key. Retrun list of (k, v) pairs.
- as_strings(alist)
Cast list elements to string.
- dict2list(adict)
Return list from dictionary, return list of (key, value) pairs.
- sort_dictionary_by_value(d, reverse=False)
Sort dictionary by value. Retrun list of (v, k) pairs.
- remove_duplicates_and_sort(data)
Remove duplicates from list and sort it.
- remove_duplicates(data)
- remove_redundancy(alist)
Remove redundancy or empty elements in given list.
- clear_fragments_redundancy(data, extend=False, same_case_func=None)
Remove nested fragments, ata format [start, end, ...].

### Ngrams (kmers) tools

	from trseeker.tools.ngrams_tools import *

- generate_ngrams(text, n=12)
Yields all ngrams of length n from given text.

- generate_window(text, n=12, step=None)
Yields all ngrams of length n from given text. 

- get_ngrams(text, m=5, n=12)
Returns m most frequent (ngram of length n, tf) tuples for given text.

- get_ngrams_freq(text, m=5, n=12)
Returns m most frequent (ngram of length n, fraction of possible ngrams) tuples for given text.

- get_ngrams_feature_set(text, m=5, n=12)
Returns a feature set {'ngram':'ngram',...}  of m most frequent ngram of length n for given text.

- get_ngram_freq_distance(ngrams_a, ngrams_b)
Returns a distance between two ngram sets where distance is a sum(min(ngram_a, ngram_b) for each common ngram). Format for ngrams_a and ngrams_b is a dictionary {ngram:n, ...}

- get_ngram_common_distance(ngrams_a, ngrams_b)
Returns a distance between two ngram sets where distance is a len(). Format for ngrams_a and ngrams_b is a dictionary {ngram:n, ...}

- get_repeatness_coefficent(length, k, kmern)
Return repeatness coefficient. From 0 (e.g. polyA) to 1 (unique sequence).

	kmern * (N-1) / (N**2)

- get_expessiveness_coefficent(kmern, k)
Return expressivenesss coefficient.

	kmern * 1. / 4^k

Update tf and df data with k-mers from given sequence.
Usage:
	
	from trseeker.tools.ngrams_tools import count_kmer_tfdf

	k = 21
	for sequence in sequences:
		tf_dict, df_dict, local_tf, local_df = count_kmer_tfdf(sequence, tf_dict, df_dict, k)


### Edit distance functions

	from trseeker.tools.edit_distance import *

- get_pos(L, d)
Return list of element (*d*) positions in given list.

- get_ed_similarity(s1, s2, monomer_mode=False, full_info=False, verbose=False)
Return edit distance between two given strings. Edit distance dynamic programming implementation. Length of second sequence does not change similarity value.
1) monomer mode: double short monomer sequence if len2/len1 < 2
2) full_info: boolean, return full information
Return: percent of ED similarity, or if full_info (distance, number of d, positions of d, S matrix)

- get_edit_distance_info(s1, s2, verbose=False, monomer_mode=False)

Return two sequences edit distance full information.    
1) s1: sequence of tandem repeat monomer
2) s2: sequence of tandem repeat monomer
3) verbose: boolean
4) monomer mode: double short monomer sequence if len2/len1 < 2.
Return: ("max_sim Nmax_sim %sim pos_list", pos_list, all_result, length seq1, length seq2)

- get_edit_distance(s1, s2)

Return ED valie for two sequences.

- get_edit_distance_row(s1, s2)

Get last row of ED matrix between two strings.

- hamming_distance(s1, s2)

Get Hamming distance: the number of corresponding symbols that differs in given strings.
    
### Working with Repbase files

	from trseeker.tools.repbase_tools import join_repbase_files

	join_repbase_files(input_folder, output_file)

- join_repbase_files(input_folder, output_file)

Function join all Repbase fasta files in one huge fasta; reformat headers for compatibility with NCBI tools.

### Working with Trace files

	from trseeker.tools.trace_tools import unclip_trace_file

	unclip_trace_file(fasta_file, clip_file, uncliped_file)

- unclip_trace_file(fasta_file, clip_file, uncliped_file)

Unclip. Remove vector flanks from Trace data.

### Function related to statistics

	from trseeker.tools.statistics import *

- get_sigma(data)

Calculate sigma, return sum(module(xi - mean).

- get_mean(data)
- get_sample_derivation(variance)
- t_test(sample_mean, dist_mean, variance, N)
- get_element_frequences(data)
- get_simple_statistics(data)

Return dictionary containing simple statistics for given list. Dictionary keys: mean, variance, sigma, sample derivation.

### Parsers

	from trseeker.tools.parsers import *

- parse_fasta_head(fa_head)
- parse_chromosome_name(head)
- trf_parse_line(line)
- trf_parse_param(line)
- trf_parse_head(line)
- get_wgs_prefix_from_ref(ref)
- get_wgs_prefix_from_head(head)

### Functions related to TRF

	from trseeker.tools.trf_tools import *

- trf_search(file_name)
- trf_search_in_dir(folder, verbose=False, file_suffix=".fa", output_folder=None)
- trf_filter_by_array_length(trf_file, output_file, cutoff)

Create output TRF file with tandem repeats with length greater than from input file. Function returns number of tandem repeats in output file.

- trf_filter_by_monomer_length(trf_file, output_file, cutoff)

Create output TRF file with tandem repeats with unit length greater than from input file. Function returns number of tandem repeats in output file.

- trf_filter_exclude_by_gi_list(trf_file, output_file, gi_list_to_exclude)

Create output TRF file with tandem repeats with GI that don't match GI_LIST. List of GI, see TRF and FA specifications, GI is first value in TRF row.

- trf_representation(trf_file, trf_output, representation)

Write TRF file tab delimited representation.Representation: numerical|index|agc_apm|with_monomer|family

- trf_write_field_n_data(trf_file, file_output, field, field_format="%s")

Write statistics data: field, N.

- trf_write_two_field_data(trf_file, file_output, field_a, field_b)

Write statistics data: field_a, field_b.

- count_trs_per_chrs(all_trf_file)

Function prints chr, all trs, 3000 trs, 10000 trs

- count_trf_subset_by_head(trf_file, head_value)

Function prints number of items with given fasta head fragment

- fix_chr_names(trf_file, temp_file_name=None, case=None)

Some fasta heads impossible to parse, so it is simpler to fix them postfactum

### Working with TRs datasets

	from trseeker.tools.trs_dataset import *

- COLORS dictionary
- get_colors(family)
- get_colors_rgb(family)
- create_mathematice_dataset_by_family(trs_dataset, path_to_mathematica_folder, min_borders, max_borders)

### Working with SRA data

	from trseeker.tools.sra_tools import *

- sra_fastaq_reader(file_name)

Iterate over fastaq data.

- sra_fasta_reader(file_name)

Iterate over fasta SRA data.

- read_fastaq_freq_to_memory(file_name)
- fastaq_to_fasta(file_name, output_file)
- write_fastaq_repeats(input_file, output_file, min_tf=1000)
- seq_to_bin(seq)
- bin_to_seq(bseq)
- write_reduced_fasta(input_file, output_reduced)

Read SRA fasta data without read repeats.

- write_ngrams(input_file, output_ngram, NGRAM_N)

Write ngrams data from SRA fasta data.

### Working with sequence patterns

	from trseeker.tools.sequence_patterns import *

Avaliable patterns:

- MARS1
- MARS2
- CENPB
- PRDB9

Functions:

- re_translate(sequence)

Translate sequence in DNA15 in reg exp.

- re_get_plus_minus(sequence)
- re_get_mutations(sequence)

Return list of sequences where one letter changed to N.

- get_double_pattern(pattern_static, pattern_dynamic)

Return list of patterns one not mutated and second mutated.

- get_mutated_pattern_twice(pattern)
- get_mutated_pattern_trice(pattern)
- pattern_search(name, sequence, pattern_function, pattern_function_params)

### Working with BLAST

	from trseeker.tools.blast_tools import *

- blastn(database, query, output)
- create_db(fasta_file, output, verbose=False, title=None)
- alias_tool(dblist, output, title)
- bl2seq(input1, input2, output)
- create_db_for_genome(file_pattern=None, chromosome_list=None, output=None, title=None)
- get_gi_list(gi_score_file, score_limit=90)
- get_all_gi_from_blast(blast_file, mode="gi")
- get_all_blast_obj_from_blast(blast_file, mode="ref")
- blastn_search_for_trs(trf_large_file, db, annotation_self_folder, temp_file, skip_by_family=None, is_huge_alpha=False)

### Working with suffix arrays

	from trseeker.tools.sa_tools import *

- fasta_to_sa_input(fasta_file, sa_input_file, index_file_name=None, file_name=None, start_id=None, increment=False)

Fasta to SA input. Text file of sequences delimeted by $ symbol.

- sa_input_to_fasta(input_file, output_file)

SA input file to fasta file.

- pickle_dictionary_for_docid_trid(sa_doc_index, doc_to_trf_file, trf_to_doc_file)

Function precompiles doc2id and id2doc pickled dictionary

- filter_sa_dataset(sa_file, output_file, min_tf, min_df)

Function filters and writes sa data with given min tf and df.

- iterate_sa_corpus(corpus)

 Yields corpus texts. Corpus is $ delimited text file.

### Working with NCBI genome information

	from trseeker.tools.ncbi_genomes import *

- load_genome_projects(subgroup)

Scrap finished WGS genome projects wtih given SubGroup name.

- load_genome_projects_exclude(notsubgroups)

Scrap finished WGS genome projects wtih given SubGroup name.

- print_add_project_data(genome_projects, pid_type, pr_type)

Print data for initiating PySatDNA projects.

### Working with graph data

	from trseeker.tools.network_tools import *

- compute_network(network_file, output_file_pattern, trf_large_index_file, g_file, weights_file, weight_to_pairs_file)

Compute network slices using graph_tool library or networkx.

- load_graph(ml_file)
- analyse_network_graph_tool(G, output_file_pattern, trf_large_index_file)
- write_classification_graph_tool(output_file, components, trid2meta)
- write_classification_neworkx(output_file, components, trid2meta)
- create_graphml(network_file, ml_file)
- load_networkx(network_file)
- init_graph_networkx(network_data, start=0, precise=1, trf_large_index_file=None)
- analyse_networkx(G, network_data, output_file_pattern, trf_large_index_file)

### Working with NCBI annotations

	from trseeker.tools.ncbi_annotation_tools import *

- get_ideogram_dict(idiogram_file, mode="NCBI")

Read ideogram file and return return dict chr -> list.

- get_gene_list(file_gene_list, color='#000000', left_padding=30, gene_group_label='C57BL/6J')

Not implemented yet.