#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 01.09.2013
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
Based on https://github.com/caballero/SeqComplex

== SeqComplex ==

This is a collection of methods to compute the composition and complexity of a DNA
sequence(s) from a Fasta file.

The SeqComplex.pm is a Perl Module with all methods, 2 examples scripts are added:
(1) compSeq.pl compute the methods in a windowed mode.
(2) profileComplexSeq.pl compute the methods using the whole sequence.

Copyright (C) 2009 by Juan Caballero [jcaballero@systemsbiology.org]

All code is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.8 or,
at your option, any later version of Perl 5 you may have available.
'''

'''
SeqComplex

=head1 SYNOPSIS

  # Example 1. Run all methods for each sequence
  use SeqComplex;
  
  my @seq = loadSeqs($file);
  my $num = 0;
  foreach my $seq (@seq) {
    $num++;
    my %results = Complex::runAllMethods( $seq, $win, $kmer );
    foreach my $m (keys %results) {
        my $res = join "\t", @{ $results{$m} };
        print "seq_$num\t$m\t$res\n";
    }
  }
  
  ## Example 2. Run a particular method 
  ## NOTE: all methods require references as inputs and returns an ARRAY.
  use SeqComplex qw(:methods)
  
  # Calculate the GC content every 1kb
  my $seq = 'ATGC' x 10000;
  my $win = 1000;
  my @gc_vals = gc(\$seq, \$win);
  
=head1 DESCRIPTION

Calculate composition and complexity of a DNA sequence.

=head2 EXPORT

runAllMethods

=head2 EXPORT_OK

:methods exports methods  = [gc at gcs ats cpg cwf ce cz ct cl cm]
:utils   exports basic utils = [log_k pot createWords countWords randSeq]
:all     exports @methods and @utils

=cut

import zlib;
k = 4 # alphabet size <A, C, G, T>
alphabet = ['A', 'G', 'C', 'T']

=head1 SUBROUTINES

=head2 Methods

=cut

=head3 runAllMethods

Function to run each method in a DNA sequence.

Call: runAllMethods ( $seq, $win, $word ) STRING, NUMBER, NUMBER

Return: %values HASH (KEY = method, VALUES = @values)

=cut

def runAllMethods {
    my $gseq   = shift @_;
    my $gwin   = shift @_; 
    $gwin    ||= length $gseq; # Use full length of the sequence
    my $gword  = shift @_;
    my @gwords = ();
    if (defined $gword) { push @gwords, $gword;    }
    else                { @gwords = (1,2,3,4,5,6); }
    
    my %gvalues = ();
    #print "Calc gc ...\n";
    @{ $gvalues{'gc' } } =  gc( \$gseq, \$gwin );
    #print "Calc at ...\n";
    @{ $gvalues{'at' } } =  at( \$gseq, \$gwin );
    #print "Calc gcs ...\n";
    @{ $gvalues{'gcs'} } = gcs( \$gseq, \$gwin );
    #print "Calc ats ...\n";
    @{ $gvalues{'ats'} } = ats( \$gseq, \$gwin );
    #print "Calc ket ...\n";
    @{ $gvalues{'ket'} } = ket( \$gseq, \$gwin );
    #print "Calc pur ...\n";
    @{ $gvalues{'pur'} } = pur( \$gseq, \$gwin );
    #print "Calc cpg ...\n";
    @{ $gvalues{'cpg'} } = cpg( \$gseq, \$gwin );
    #print "Calc ce ...\n";
    @{ $gvalues{'ce' } } =  ce( \$gseq, \$gwin );
    #print "Calc cz ...\n";
    @{ $gvalues{'cz' } } =  cz( \$gseq, \$gwin );
    #print "Calc cwf ...\n";
    @{ $gvalues{'cwf'} } = cwf( \$gseq, \$gwin );
    foreach my $ws (@gwords) {
    #   print "Calc ct$ws ...\n";
        @{ $gvalues{"ct$ws" } } =  ct( \$gseq, \$gwin, \$ws );
    #   print "Calc cl$ws ...\n";
        @{ $gvalues{"cl$ws" } } =  cl( \$gseq, \$gwin, \$ws );
    #   print "Calc cm$ws ...\n";
        @{ $gvalues{"cm$ws" } } =  cm( \$gseq, \$gwin, \$ws );
    }
    return %gvalues;
}

=head3 gc

Function to calculate the GC content.

Call: gc( \$seq, \$win ) STRING, NUMBER

Return: @values ARRAY

=cut

sub gc {
    my $seq    = shift;
    my $win    = shift;
    my $len    = length $$seq;
    my @values = ();
    for (my $p = 0; $p <= ($len - $$win); $p += $$win) {
        my $str = substr ($$seq, $p, $$win);
        my %elm = countWords($str, 1);
        my $r   = 0;
        my $tot = $elm{'C'} + $elm{'G'} + $elm{'T'} + $elm{'A'};
        $r      = ($elm{'C'} + $elm{'G'}) / $tot if($tot > 1);
        push @values, $r; 
    }
    return @values;
}

=head3 at

Function to calculate the AT content.

Call: at( \$seq, \$win ) STRING, NUMBER

Return: @values ARRAY

=cut

sub at {
    my $seq    = shift;
    my $win    = shift;
    my $len    = length $$seq;
    my @values = ();
    for (my $p = 0; $p <= ($len - $$win); $p += $$win) {
        my $str = substr ($$seq, $p, $$win);
        my %elm = countWords($str, 1);
        my $r   = 0;
        my $tot = $elm{'C'} + $elm{'G'} + $elm{'T'} + $elm{'A'};
        $r      = ($elm{'A'} + $elm{'T'}) / $tot if($tot > 1);
        push  @values, $r;
    }
    return @values;
}

=head3 gcs

Function to calculate the GC skew content.

Call: gcs( \$seq, \$win ) STRING, NUMBER

Return: @values ARRAY

=cut

sub gcs {
    my $seq    = shift;
    my $win    = shift;
    my $len    = length $$seq;
    my @values = ();
    for (my $p = 0; $p <= ($len - $$win); $p += $$win) {
        my $str = substr ($$seq, $p, $$win);
        my %elm = countWords($str, 1);
        my $l   = $elm{'G'} + $elm{'C'};
        my $r   = 0;
        $r      = abs($elm{'G'} - $elm{'C'}) / $l if ($l > 1);
        push  @values, $r; 
    }
    return @values;
}

=head3 cpg

Function to calculate the CpG content.

Call: cpg( \$seq, \$win ) STRING, NUMBER

Return: @values ARRAY

=cut

sub cpg {
    my $seq    = shift;
    my $win    = shift;
    my $len    = length $$seq;
    my @values = ();
    for (my $p = 0; $p <= ($len - $$win); $p += $$win) {
        my $str = substr ($$seq, $p, $$win);
        my %elm = countWords($str, 1);
        my $c   = $str =~ tr/CG/CG/;
        my $l   = $elm{'G'} * $elm{'C'};
        my $r   = 0;
        $r      = $c / $l if ($l > 1);
        push @values, $r; 
    }
    return @values;
}

=head3 ats

Function to calculate the AT skew content.

Call: ats( \$seq, \$win ) STRING, NUMBER

Return: @values ARRAY

=cut

sub ats {
    my $seq    = shift;
    my $win    = shift;
    my $len    = length $$seq;
    my @values = ();
    for (my $p = 0; $p <= ($len - $$win); $p += $$win) {
        my $str = substr ($$seq, $p, $$win);
        my %elm = countWords($str, 1);
        my $l   = $elm{'A'} + $elm{'T'};
        my $r   = 0;
        $r      = abs($elm{'A'} - $elm{'T'}) / $l if ($l > 1);
        push  @values, $r;
    }
    return @values;
}

=head3 ket

Function to calculate the Keto skew content.

Call: ket( \$seq, \$win ) STRING, NUMBER

Return: @values ARRAY

=cut

sub ket {
    my $seq    = shift;
    my $win    = shift;
    my $len    = length $$seq;
    my @values = ();
    for (my $p = 0; $p <= ($len - $$win); $p += $$win) {
        my $str = substr ($$seq, $p, $$win);
        my %elm = countWords($str, 1);
        my $r   = 0;
        my $tot = $elm{'G'} + $elm{'T'} + $elm{'A'} + $elm{'C'};
        $r      = abs($elm{'G'} + $elm{'T'} - $elm{'A'} - $elm{'C'}) / $tot if($tot > 1);
        push @values, $r; 
    }
    return @values;
}

=head3 pur

Function to calculate the Purine skew content.

Call: pur( \$seq, \$win ) STRING, NUMBER

Return: @values ARRAY

=cut

sub pur {
    my $seq    = shift;
    my $win    = shift;
    my $len    = length $$seq;
    my @values = ();
    for (my $p = 0; $p <= ($len - $$win); $p += $$win) {
        my $str = substr ($$seq, $p, $$win);
        my %elm = countWords($str, 1);
        my $r   = 0;
        my $tot = $elm{'G'} + $elm{'T'} + $elm{'A'} + $elm{'C'};
        $r      = abs($elm{'G'} - $elm{'T'} + $elm{'A'} - $elm{'C'}) / $tot if($tot > 1);
        push @values, $r; 
    }
    return @values;
}

=head3 cwf

Function to calculate the Complexity by Wootton & Federhen values.

Call: cwf( \$seq, \$win, ) STRING, NUMBER

Return: @values ARRAY

=cut

sub cwf {
    my $seq    = shift;
    my $win    = shift;
    my $len    = length $$seq;
    my @values = ();
    my $up = 0;
    for (my $i = 1; $i <= $$win; $i++) { $up += log_k($k, $$win); }
    for (my $p = 0; $p <= ($len - $$win); $p += $$win) {
        my $str = substr ($$seq, $p, $$win);
        my %elm = countWords($str, 1);
        my $r   = 0;
        my $dw  = 0;
        my $tot = $elm{'G'} + $elm{'T'} + $elm{'A'} + $elm{'C'};

        foreach my $b (keys %elm) {
            next unless ($elm{$b} > 0);
            $dw += log_k($k, $elm{$b}); 
        }
        $r = ($up - $dw) / $tot if($tot > 1);;
        push @values, $r; 
    }
    return @values;
}

=head3 ce

Function to calculate the Complexity Entropy values.

Call: ce( \$seq, \$win ) STRING, NUMBER

Return: @values ARRAY

=cut

sub ce {
    my $seq    = shift;
    my $win    = shift;
    my $len    = length $$seq;
    my @values = ();
    for (my $p = 0; $p <= ($len - $$win); $p += $$win) {
        my $str = substr ($$seq, $p, $$win);
        my %elm = countWords($str, 1);
        my $ce  = 0;        
        my $tot = $elm{'G'} + $elm{'T'} + $elm{'A'} + $elm{'C'};

        foreach my $b (keys %elm) {
            next unless ($elm{$b} > 0);
            my $r = 0; 
            $r    = $elm{$b} / $tot if($tot > 1); 
            $ce  -= $r * log_k($k, $r); 
        }
        push @values, $ce;
    }
    return @values
}

=head3 cm

Function to calculate the Complexity in Markov model values.

Call: cm( \$seq, \$win, \$word ) STRING, NUMBER, NUMBER

Return: @values ARRAY

=cut

sub cm {
    my $seq    = shift;
    my $win    = shift;
    my $word   = shift;
    my $len    = length $$seq;
    my @values = ();
    my $m = pot($k, $$word);
    for (my $p = 0; $p <= ($len - $$win); $p += $$win) {
        my $str = substr ($$seq, $p, $$win);
        my %elm = countWords($str, $$word);
        my $cm  = 0;
        my $dw  = $$win - $$word - 1;
        foreach my $b (keys %elm) {
            next unless ($elm{$b} > 0);
            my $r = $elm{$b} / $dw; 
            $cm  -= $r * log_k($m, $r); 
        }
        push @values, $cm;
    }
    return @values;
}

=head3 cl

Function to calculate the Complexity Linguistic values.

Call: cl( \$seq, \$win, \$word ) STRING, NUMBER, NUMBER

Return: @values ARRAY

=cut

sub cl {
    my $seq    = shift;
    my $win    = shift;
    my $word   = shift;
    my $len    = length $$seq;
    my @values = ();
    for (my $p = 0; $p <= ($len - $$win); $p += $$win) {
        my $str     = substr ($$seq, $p, $$win);
        my $sum_vl  = 0;
        my $sum_vm  = 0;
        for (my $l = 1; $l <= $$word; $l++) {
            my $vl  = 0;
            my $vm  = 0;
            my $pot = pot($k, $l); 
            if ($pot < ($$win - $l + 1)) { $vm = $pot;           }
            else                         { $vm = $$win - $l + 1; }
            $sum_vm += $vm;
            my %elm = countWords($str, $l);
            foreach my $b (keys %elm) {
                next unless ($elm{$b} > 0);
                $vl++; 
            }
            $sum_vl += $vl;
        }
        my $r = 0;
        $r    = $sum_vl / $sum_vm if ($sum_vm > 0);
        push @values, $r;
    }
    return @values;
}

=head3 ct

Function to calculate the Complexity by Trifonov values.

Call: ct( \$seq, \$win, \$word ) STRING, NUMBER, NUMBER

Return: @values ARRAY

=cut

sub ct {
    my $seq    = shift;
    my $win    = shift;
    my $word   = shift;
    my $len    = length $$seq;
    my @values = ();
    for (my $p = 0; $p <= ($len - $$win); $p += $$win) {
        my $ct  = 1;
        my $str = substr ($$seq, $p, $$win);
        for (my $l = 1; $l <= $$word; $l++) {
            my $vl  = 0;
            my $vm  = 0;
            my $pot = pot($k, $l);
            if ($pot < ( $$win - $l + 1 )) { $vm = $pot;           }
            else                           { $vm = $$win - $l + 1; }
            my %elm = countWords($str, $l);
            foreach my $b (keys %elm) {
                next unless ($elm{$b} > 0);
                $vl++; 
            }
            $ct *= $vl / $vm;
        }
        push @values, $ct;
    }
    return @values;
}

=head3 cz

Function to calculate the Compression factor.

Call: cz( \$seq, \$win ) STRING, NUMBER

Return: @values ARRAY

=cut

sub cz {
    my $seq    = shift;
    my $win    = shift;
    my $len    = length $$seq;
    my @values = ();
    for (my $p = 0; $p <= ($len - $$win); $p += $$win) {
        my $str = substr ($$seq, $p, $$win);
        my $r   = 0;
        my $z   = Compress::Zlib::memGzip($str);
        my $bz  = length $z;
        $r      = $bz / $$win; 
        push @values, $r;
    }
    return @values;
}

=head3 clz

Function to calculate the CLZ values.

Call: clz( \$seq, \$win ) STRING, NUMBER

Return: @values ARRAY

* Not implemented yet!

=cut

sub clz {
    print "Sorry, not implemented yet!\n";
    exit;
}

=head2 UTILITIES

=head3 pot

Function for calculate the exponential of a number.

Call: pot( $num, $exp ) NUMBER, NUMBER

Return: $res NUMBER

=cut

sub pot {
    my $num = shift @_;
    my $exp = shift @_;
    if ($num == 4) {
        if    ($exp == 1) { return       4; }
        elsif ($exp == 2) { return      16; }
        elsif ($exp == 3) { return      64; }
        elsif ($exp == 4) { return     256; }
        elsif ($exp == 5) { return    1024; }
        elsif ($exp == 6) { return    4096; }
        elsif ($exp == 7) { return   16384; }
        elsif ($exp == 8) { return   65536; }
        elsif ($exp == 9) { return  262144; }
        elsif ($exp ==10) { return 1048576; }
    }
    my $res = $num;
    for (my $i = 2; $i <= $exp; $i++) { $res *= $num; }
    return $res;
}

=head3 log_k

Function for calculate the logarithm of any base.

Call: log_k( $base, $num ) NUMBER, NUMBER

Return: $res NUMBER

=cut

sub log_k {
    my $base = shift @_;
    my $num  = shift @_;
    my $res  = 0;
    if ($num > 0) { 
        $res = log $num / log $base; 
    }
    else {
        die "Cannot calculate log_k($base, $num)\n";
    }
    return $res;
}

=head3 countWords

Function for count words in a sequence.

Call: countWords( $seq, $word ) STRING, NUMBER

Return: %results HASH (KEYS are the elements)

=cut

sub countWords {
    my $seq   = shift @_;
    my $word  = shift @_;
    my $len   = length $seq; 
    my %count = ();
    
    # Init state when word == 1
    foreach my $b (@alphabet) {
        $count{$b} = 0;
    }
    
    for (my $i = 0; $i <= ($len - $word); $i++) {
        my $elem = substr ($seq, $i, $word);
        next if($elem =~ /[^ACGT]/);
        $count{$elem}++;
    }
    return %count;
}
'''
import sys, os
try:
    import PyExp
except:
    sys.path.append("/Users/ad3002/Dropbox/workspace")
    sys.path.append("/root/Dropbox/workspace")
    sys.path.append("c:/Users/Master/Dropbox/workspace")
import zlib

def get_zlib_complexity(s):
    ''' Get simple complexity by zlib compression ratio.
    '''
    return float(len(zlib.compress(s)))/len(s)

def createWords():
    '''Function to create a list of possible words with a specific length
        Call: createWords( $kmer, @alphabet ) NUMBER, ARRAY

        Return: @old ARRAY
    '''
    '''
    my $kmer = shift @_; $kmer--;
    my @old = @_;
    my @new = ();
    if ($kmer < 1) {
        return @old;
    }
    else {
        foreach my $e (@old) {
            foreach my $n (@alphabet) {
                push @new, "$e$n"; # add new element
            }
        }
        createWords($kmer, @new); # recursion call
    }
    '''
    pass

def get_rand_seq():
    ''' Function to create a random DNA sequence'''
    pass

if __name__ == '__main__':
    

    from trseeker.seqio.fasta_file import sc_iter_fasta
    from collections import defaultdict
    from trseeker.tools.ngrams_tools import process_list_to_kmer_index

    S = []
    for i, seq_obj in enumerate(sc_iter_fasta("/root/Dropbox/rna.fa")):
        if i > 4000:
            break
        S.append(seq_obj.sequence)
    d = process_list_to_kmer_index(S, 23)
    total = 0
    gH = defaultdict(int)
    for x in d:
        gH[get_zlib_complexity(x[0])] += int(x[2])
        total += float(x[2])
    for k in gH:
        gH[k] = int(1000*gH[k]/total)

    S = []
    for i, seq_obj in enumerate(sc_iter_fasta("/home/repbase/repbase.fa")):
        if i > 4000:
            break
        S.append(seq_obj.sequence)
    d = process_list_to_kmer_index(S, 23)
    total = 0
    rH = defaultdict(int)
    for x in d:
        rH[get_zlib_complexity(x[0])] += int(x[2])
        total += float(x[2])
    for k in rH:
        rH[k] = int(1000*rH[k]/total)



    fh = open("/home/animals_wgs/botryllus_schlosseri_wgs/trf_large_kmers.index")
    d = fh.readlines()
    d = [x.split("\t") for x in d]
    sH = defaultdict(int)
    cH = defaultdict(int)
    total = 0
    for x in d:
        sH[get_zlib_complexity(x[0])] += int(x[2])
        total += float(x[2])
    for k in sH:
        sH[k] = int(1000*sH[k]/total)

    fh = open("/home/animals_wgs/botryllus_schlosseri_wgs/trf_micro_kmers.index")
    d = fh.readlines()
    d = [x.split("\t") for x in d]
    total = 0
    for x in d:
        cH[get_zlib_complexity(x[0])] += int(x[2])
        total += float(x[2])
    for k in cH:
        cH[k] = int(1000*cH[k]/total)

    with open("/root/Dropbox/trf_large_kmers.hist", "w") as fh:
        for k in sH:
            fh.write("%s\t%s\t%s\t%s\t%s\n" % (k, sH[k], cH[k], gH[k], rH[k]))
        for k in cH:
            if k in sH:
                continue
            fh.write("%s\t%s\t%s\t%s\t%s\n" % (k, sH[k], cH[k], gH[k], rH[k]))
        for k in gH:
            if k in sH or k in cH:
                continue
            fh.write("%s\t%s\t%s\t%s\t%s\n" % (k, sH[k], cH[k], gH[k], rH[k]))
        for k in rH:
            if k in sH or k in cH or k in gH:
                continue
            fh.write("%s\t%s\t%s\t%s\t%s\n" % (k, sH[k], cH[k], gH[k], rH[k]))