#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2007-2009 Aleksey Komissarov ( ad3002@gmail.com )
# All rights reserved.
#
# This software is licensed as described in the file COPYING, which
# you should have received as part of this distribution.
"""
Function for simple sequence analysis.

- get_revcomp(sequence)
- get_gc(sequence)
- check_gapped(sequence)
- clear_sequence(sequence)
- restriction(sequence, pattern, end="")
- restriction_fragments_n(sequence, pattern)
- get_subseq(seq, start, end)
- check_synonims(seq_a, seq_b)
- random_mutation(seq, n, alphabet="actgn ")
    
TODO: Not implemented:

- get_bin_sequence(sequence)
    
"""
import re
import random
from trseeker.tools.other_tools import sort_dictionary_by_value
from collections import defaultdict

def get_shifts_variants(sequence):
    '''
    '''
    shifts = set()
    for i in xrange(len(sequence)):
        shifts.add(sequence[i:]+sequence[:i])
    return list(shifts)
    
def get_revcomp(sequence):
    '''Return complementary sequence.

    >>> complementary('AT CG')
    'CGAT'

    '''
    __complementary = dict(zip('ATCGNatcgn[]', 'TAGCNtagcn]['))
    c = __complementary
    return ''.join(c.get(nucleotide, '') for nucleotide in reversed(sequence))

def fix_strand(sequence):
    ''' Return normalized sequence with following rules:
    1) if T > A then return reverse complement
    2) if A == T and G > C  then return reverse complement
    3) else return sequence
    '''
    A = sequence.count('a')
    T = sequence.count('t')
    if T > A:
        return get_revcomp(sequence)
    if T == A:
        C = sequence.count('c')
        G = sequence.count('g')
        if G > C:
            return get_revcomp(sequence)
    return sequence

def get_gc(sequence):
    ''' Count GC content.
    '''
    length = len(sequence)
    if not length:
        return 0
    sequence = clear_sequence(sequence)
    count_c = sequence.count('c')
    count_g = sequence.count('g')
    gc = float(count_c + count_g) / float(length)
    return float(gc)
    
def check_gapped(sequence):
    """ Check for N in sequence. Return n(with N) or w(whole).
    """
    w_regexp = re.compile('n|N')
    regexp_obj = w_regexp.search(sequence)
    if (regexp_obj):
        return True
    else:
        return False
    
def get_subseq(seq, start, end):
    ''' Return subsequence.
    '''
    return seq[start-1 : end]
    
def clear_sequence(sequence):
    """ Clear sequence (ACTGN alphabet):
    
    - lower case
    - \s -> ""
    - [^actgn] -> ""
    """
    sequence = sequence.lower().strip()
    sequence = re.sub("\s+", "", sequence)
    return re.sub("[^actgn]", "", sequence)

def restriction(sequence, pattern, end=None):
    ''' Return list of sequences splited by pattern with added end.
    
    - patterm: reg_exp
    '''
    fragments = [ fragment for fragment in re.split(pattern, sequence) ]
    if len(fragments) == 1 and \
        fragments[0] == sequence:
            return fragments
    n = len(fragments)
    if not end:
        end = ""
    for i in range(0,n):
        if i == 0:
            fragments[i] += end
        elif i == n-1:
            fragments[i] = end + fragments[i]
        else:
            fragments[i] = end + fragments[i] + end
    return fragments
            
def restriction_fragments_n(sequence, pattern):
    '''Return number of fragments after restriction of sequence with given pattern.
    
    - patterm: reg_exp
    '''
    fragments = [ fragment for fragment in re.split(pattern, sequence) ]
    n = len(fragments)
    m = [fragment for fragment in fragments if fragment == ""]    
    return n - len(m)

def check_cyclic_repeats(seq_a, seq_b):
    ''' Check tandem repeat synonims  between two sequence with same length.

    >>> check_synonims("acc", "accacc")
    True 
    
    @param seq_a: sequence with length n
    '''
    if len(seq_a) != len(seq_b):
        print "Need sequences with same length"
        return False
    if seq_a in seq_b*2 or seq_b in seq_a*2:
        return True
    else:
        return False

def random_mutation(seq, n, alphabet="actgn +"):
    ''' Return sequence with n mutations.
    
    - alphabet: letters and ' ' - del, '+' - inserts
    '''
    N = len(seq)-1
    seq = [ x[1] for x in enumerate(seq) ]
    for mut in range(0, n):
        i = random.randint(0, N)
        j = random.randint(0, len(alphabet)-1)
        seq[i] = alphabet[j]
    mutated = "".join(seq)
    mutated.replace(" ","")
    mutated.replace("+",alphabet[random.randint(0, len(alphabet)-1)])
    return mutated

def get_consensus(strs):
    ''' Return consensus string for given list of strings.'''
    n = len(strs[0])
    for s in strs:
        assert len(s) == n
    result = {}
    result['A'] = [0]*n
    result['C'] = [0]*n
    result['T'] = [0]*n
    result['G'] = [0]*n
    for i in xrange(n):
        for s in strs:
            result[s[i]][i] += 1
    return result

def remove_consensus_redundancy(trf_objs):
    ''' Take a minimal sequence from lexicographically sorted rotations of sequence and its reverse complement
        Example: ACT, underlined - reverse complement seqeunces
        ACT, AGT, CTA, GTA, TAC, TAG
        Find all possible multimers, e.g. replace GTAGTAGTA consensus sequence with ACT
        Return:
        1) list of sorted TRs
        2) list of (df, consensus) pairs
    '''
    # sort by length
    consensuses = [x.trf_consensus for x in trf_objs]
    consensuses = list(set(consensuses))
    consensuses.sort(key=lambda x: len(x))
    length2consensuses = defaultdict(list)
    result_rules = {}
    for i, monomer in enumerate(consensuses):
        n = len(monomer)
        length2consensuses[n].append(i)
    for i, monomer in enumerate(consensuses):
        if monomer in result_rules:
            continue
        base = len(monomer)
        n = base
        # maximal consensus length from TRF is 2000 bp
        variants = set(get_shifts_variants(monomer) + get_shifts_variants(get_revcomp(monomer)))
        lex_consensus = min(variants)        
        while n <= 2000:
            for k in length2consensuses[n]:
                monomer_b = consensuses[k]
                if monomer_b in result_rules:
                    continue
                s = n/base
                v = set()
                for p in xrange(s):
                    v.add(monomer_b[p*base:(p+1)*base])            
                if len(v) > 1:
                    continue
                item = v.pop()
                if item in variants:
                    if consensuses[k] != lex_consensus:
                        print i, base, consensuses[k], "->", lex_consensus, "\r",
                    result_rules[consensuses[k]] = lex_consensus
                    
            n += base

    variants2df = defaultdict(int)
    for trf_obj in trf_objs:
        variants2df[result_rules[trf_obj.trf_consensus]] += 1
        trf_obj.trf_consensus = result_rules[trf_obj.trf_consensus]
        
    print
    print "Sort families by df..."
    variants2df = sort_dictionary_by_value(variants2df, reverse=True)
    return trf_objs, variants2df
