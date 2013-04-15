#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 23.02.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
Functions related to ngram (k-mer).
'''
from collections import defaultdict
from trseeker.tools.sequence_tools import fix_strand, get_revcomp
from trseeker.tools.edit_distance import hamming_distance
from trseeker.seqio.tab_file import sc_iter_tab_file
from trseeker.models.trf_model import TRModel

def generate_ngrams(text, n=12):
    ''' Yields all ngrams of length k from given text. 
    
    - n: ngram length
    '''
    for i in xrange(0, len(text) - n + 1):
        yield i, text[i:i + n]

def generate_window(text, n=12, step=None):
    ''' Yields all ngrams of length n from given text. 
    - n: ngram length
    - step: window step
    '''
    if not step:
        step = n / 2
    for i in xrange(0, len(text) - n + 1, step):
        yield i, text[i:i + n]

def get_ngrams(text, m=5, n=12):
    ''' Returns m most frequent (ngram of length n, tf) tuples for given text.
    
    - m: number of returned ngrams
    - n: ngram length
    '''

    ngrams = {}
    for pos, ngram in generate_ngrams(text, n=n):
        ngrams.setdefault(ngram, 0)
        ngrams[ngram] += 1
    ngrams = [(key, value) for key, value in ngrams.items()]
    ngrams.sort(reverse=True, key=lambda x: x[1])
    return ngrams[:m]

def get_ngrams_freq(text, m=500000, n=12):
    ''' Returns m most frequent (ngram of length n, fraction of possible ngrams) tuples for given text.
    
    - m: number of returned ngrams
    - n: ngram length
    '''

    ngrams = defaultdict(int)
    for pos, ngram in generate_ngrams(text, n=n):
        ngrams[ngram] += 1
    text_length = float(len(text) - n + 1)
    ngrams = [(key, value, value / text_length) for key, value in ngrams.items()]
    ngrams.sort(reverse=True, key=lambda x: x[1])
    return ngrams[:m]

def get_ngrams_feature_set(text, m=5, n=12):
    ''' Returns a feature set {'ngram':'ngram',...}  of m most frequent ngram of length n for given text.
    
    - m: number of returned ngrams
    - n: ngram length
    '''

    ngrams = {}
    for pos, ngram in generate_ngrams(text, n=n):
        ngrams.setdefault(ngram, 0)
        ngrams[ngram] += 1
    result = [(key, value) for key, value in ngrams.items()]
    ngrams.sort(reverse=True, key=lambda x: x[1])
    data = {}
    for key, value in enumerate(result[:m]):
        data[key] = key
    return data

def get_ngram_freq_distance(ngrams_a, ngrams_b):
    ''' Returns a distance between two ngram sets where distance is a sum(min(ngram_a, ngram_b) for each common ngram)
    
    - ngrams_a: dictionary {ngram:n, ...}
    - ngrams_b: dictionary {ngram:n, ...}
    '''
    distance = 0
    common_ngrams = [ngram for ngram, fr in ngrams_a.items() if ngram in ngrams_b]
    for ngram in common_ngrams:
        distance += min(ngrams_a[ngram], ngrams_b[ngram])
    return distance

def get_ngram_common_distance(ngrams_a, ngrams_b):
    ''' Returns a distance between two ngram sets where distance is a len()
    
    - ngrams_a: dictionary {ngram:n, ...}
    - ngrams_b: dictionary {ngram:n, ...}
    '''
    common_ngrams = [ngram for ngram, fr in ngrams_a.items() if ngram in ngrams_b]
    return len(common_ngrams)

def get_repeatness_coefficent(length, k, kmern):
    ''' Return repeatness coefficient. From 0 (e.g. polyA) to 1 (unique sequence).
    '''
    N = length - k + 1.
    return kmern * (N-1) / (N**2)

def get_expessiveness_coefficent(kmern, k):
    ''' Return expressivenesss coefficient.
    '''
    return kmern * 1. / 4^k

def count_kmer_tfdf(sequence, tf_dict, df_dict, k):
    ''' Update tf and df data with k-mers from given sequence.
    '''
    seen = set()
    local_tf = defaultdict(int)
    local_df = defaultdict(int)
    for (ngram, tf, nf) in get_ngrams_freq(sequence, m=100000, n=k):
        if 'n' in ngram:
            continue
        seen.add(ngram)
        tf_dict[ngram] += tf
        local_tf[ngram] += tf
    for ngram in seen:
        df_dict[ngram] += 1
        local_df[ngram] += 1
    return tf_dict, df_dict, local_tf, local_df

def get_kmer_tf_df_for_data(data, k, docids=False):
    '''
    '''
    df = defaultdict(int)
    tf = defaultdict(int)
    kmer2ids = defaultdict(list)
    kmer2freq = defaultdict(list)
    verbose = False
    N = len(data)
    if N>100:
        verbose = True
    for i, sequence in enumerate(data):
        if verbose:
            print "Process td/df: ", i, N, "\r", 
        tf, df, local_tf, local_df = count_kmer_tfdf(sequence, tf, df, k)
        if docids:
            for key in local_tf:
                kmer2ids[key].append(i)
                kmer2freq[key].append(local_tf[key])
    print
    if docids:
        return (tf, df, kmer2ids, kmer2freq)
    return (tf, df)
    
def get_df_stats_for_list(data, k, kmer2df):
    ''' Compute max df, number and procent of sequence with given ngram.
    Return (maxdf, nmaxdf, pmaxdf)
    '''
    df = defaultdict(int)
    tf = defaultdict(int)
    n = len(data)
    ngram_seqs = []
    for sequence in data:
        tf, df, local_tf, local_df = count_kmer_tfdf(sequence, tf, df, k)
    result = [(v,k) for (k,v) in df.items()]
    result.sort()
    maxdf = result[-1][0]
    ngram_seqs = [(k, tf[k], kmer2df[k]) for v,k in result if v == maxdf]
    ngram_seqs.sort(key=lambda x: x[1], reverse=True)
    nmaxdf = len(ngram_seqs)
    pmaxdf = round(float(maxdf)/n, 3)
    ngram_seqs = [":".join((k,str(f), str(d))) for k,f,d in ngram_seqs[:10]]
    return (maxdf, nmaxdf, pmaxdf, ngram_seqs)

def process_list_to_kmer_index(data, k, docids=True, cutoff=None):
    ''' Get list of string.
    Return list of (kmer, revkmer, tf, df, docids)
    '''
    if docids:
        (tf_dict, df_dict, doc_data, freq_data) = get_kmer_tf_df_for_data(data, k, docids=docids)
    else:
        (tf_dict, df_dict) = get_kmer_tf_df_for_data(data, k, docids=docids)
    result = []
    seen = set()
    for key in df_dict:
        if key in seen:
            continue
        revkey = get_revcomp(key)
        if revkey in seen:
            continue
        if revkey in df_dict:
            df = df_dict[key] + df_dict[revkey]
            tf = tf_dict[key] + tf_dict[revkey]
            if docids:
                ids = doc_data[key] + doc_data[revkey]
                freqs = freq_data[key] + freq_data[revkey]
        else:
            df = df_dict[key]
            tf = tf_dict[key]
            if docids:
                ids = doc_data[key]
                freqs = freq_data[key]
        if docids:
            result.append((key, revkey, tf, df, ids, freqs))
        else:
            result.append((key, revkey, tf, df, "", ""))
        seen.add(key)
        seen.add(revkey)
    if cutoff:
        result = [x for x in result if x[-3]>cutoff]
    result.sort(key=lambda x: x[-3], reverse=True)
    return result

def compute_kmer_index_for_fasta_file(file_name, index_file, k=23):
    """ 
    """
    data = []
    print "Read arrays..."
    for i, seq_obj in enumerate(sc_iter_fasta(file_name)):
        data.append(seq_obj.sequence)
    print "Readed %s arrays." % i
    print "Compute k-mers..."
    result = process_list_to_kmer_index(data, k, docids=False)
    print "Save index..."
    with open(index_file, "w") as fh:
        for item in result:
            s  = "%s\n" % "\t".join(map(str, item))
            fh.write(s)
    return result

def compute_kmer_index_for_trf_file(file_name, index_file, k=23):
    """ 
    """
    data = []
    print "Read arrays..."
    for i, trf_obj in enumerate(sc_iter_tab_file(file_name, TRModel)):
        data.append(trf_obj.trf_array)
    print "Readed %s arrays." % i
    print "Compute k-mers..."
    result = process_list_to_kmer_index(data, k, docids=False)
    print "Save index..."
    with open(index_file, "w") as fh:
        for item in result:
            s  = "%s\n" % "\t".join(map(str, item))
            fh.write(s)
    return result

def get_sequence_kmer_coverage(sequence, kmers, k):
    '''
    '''
    n = len(sequence)
    match = 0.
    mismatch = 0.
    variability = set()
    for i, kmer in generate_ngrams(sequence, n=k):
        if kmer in kmers:
            variability.add(kmer)
            match += 1
        else:
            mismatch += 1
    return match/(n-k+1), variability

