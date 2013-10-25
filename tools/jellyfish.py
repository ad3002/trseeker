'''
Jellyfish python wrapper
'''
import os
from trseeker.settings import load_settings
from PyExp import sc_iter_filepath_folder
import subprocess
from trseeker.tools.seqfile import sort_file_by_int_field
from trseeker.tools.ngrams_tools import process_list_to_kmer_index
from trseeker.tools.sequence_tools import get_revcomp

settings = load_settings()
location = settings["blast_settings"]["jellyfish_location"]

def count_kmers(input_file, ouput_prefix, k, mintf=None):
    '''Count kmers with Jellyfish.

    Args:
        input_fasta: An input fasta or fastq file name.
        ouput_prefix: An output path with file prefix.
        k: A k-mer length.
        mintf: Count only k-mers with frequency greater than mintf.

    Returns:
        None.
    '''
    params = {
        "location": location,
        "input_fasta": input_file,
        "k": k,
        "hash_size": settings["jellyfish_settings"]["hash_size"],
        "hash_bits": settings["jellyfish_settings"]["hash_bits"],
        "threads": settings["jellyfish_settings"]["threads"],
        "both_strands": settings["jellyfish_settings"]["both_strands"],
        "ouput_prefix": ouput_prefix,
        "mintf": "",
    }
    if mintf:
        params["mintf"] = "--lower-count=%s" % mintf
    command = "%(location)s count %(mintf)s -m %(k)s -o %(ouput_prefix)s -c %(hash_bits)s -s %(hash_size)s %(both_strands)s -t %(threads)s %(input_fasta)s" % params
    print "Execute:", command
    return os.system(command)

def merge_kmers(folder, ouput_prefix, ouput_file):
    '''Merge Jellyfish count output to one file.
    '''
    if not ouput_file.endswith(".jf"):
        ouput_file += ".jf"
    params = {
        "location": location,
        "ouput_file": ouput_file,
        "ouput_prefix": ouput_prefix,
    }
    file_count = 0
    for file_name in sc_iter_filepath_folder(folder, mask="."):
        if ouput_prefix in file_name:
        # head, tail = os.path.split(ouput_file)
        # if tail in file_name:
            file_count += 1
    assert file_count>0
    if file_count == 1:
        command = "cp %(ouput_prefix)s_0 %(ouput_file)s" % params
    else:
        command = "%(location)s merge -o %(ouput_file)s %(ouput_prefix)s\_*" % params
    print "Execute:", command
    os.system(command)

    command = "rm %(ouput_prefix)s_*" % params
    print command
    os.system(command)

def stats_kmers(db_file, stats_file):
    '''Compute statistics.
    '''
    params = {
        "location": location,
        "db_file": db_file,
        "stats_file": stats_file,
    }
    command = "%(location)s stats --verbose -o %(stats_file)s %(db_file)s" % params
    print "Execute:", command
    os.system(command)

def histo_kmers(db_file, histo_file):
    '''Compute frequence histogram.
    '''
    params = {
        "location": location,
        "db_file": db_file,
        "histo_file": histo_file,
    }
    command = "%(location)s histo --verbose -o %(histo_file)s %(db_file)s" % params
    print "Execute:", command
    os.system(command)

def dump_kmers(db_file, kmers_file, dumpmintf):
    '''Dump k-mers database to tab-delimited file.
    '''
    params = {
        "location": location,
        "db_file": db_file,
        "kmers_file": kmers_file,
        "dumpmintf": dumpmintf,
    }
    command = "%(location)s dump --column -L %(dumpmintf)s --tab -o %(kmers_file)s %(db_file)s" % params
    print "Execute:", command
    os.system(command)
    sort_file_by_int_field(kmers_file, 1)

def query_kmers(db_file, query_hashes, both_strands=True, verbose=True):
    '''Query jellyfish database.
    @return: dictionary hash to tf
    '''
    params = {
        "location": location,
        "db_file": db_file,
        "query_hashes": query_hashes,
        "both_strands": "",
    }
    if both_strands:
        params["both_strands"] = "-C"
    command = "%(location)s query %(both_strands)s %(db_file)s" % params
    if verbose:
        print command
    final_result = {}
    n = len(query_hashes)
    for k in xrange(0,n,100):
        pp = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE, universal_newlines=True)
        for i, query in enumerate(query_hashes[k:k+100]):
            pp.stdin.write(query+" ")
        data = pp.communicate()
        error = data[1]
        if "Can't open file" in error:
            return None
        for item in data[0].strip().split("\n"):
            if item:
                if verbose:
                    print item
                key, value = item.strip().split()
                final_result[key] = value
            else:
                if verbose:
                    print data
                final_result[-1] = data[1]
    return final_result

def get_kmer_db_and_fasta(folder, input_file, kmers_file, k=23, mintf=None):
    '''Count kmer, merge them, and save to tab-delimited kmers_file.
    '''
    count_kmers(input_file, input_file, k, mintf=mintf)
    merge_kmers(folder, input_file, input_file)
    db_file = "%s.jf" % input_file
    dump_kmers(db_file, kmers_file)

def query_and_write_coverage_histogram(db_file, query_sequence, output_file, k=23):
    '''Save coverage histogram into output_file for given query_sequence.
    '''
    index = process_list_to_kmer_index([query_sequence], k, docids=False, cutoff=None)
    query_hashes = [x[0] for x in index]
    data =  query_kmers(db_file, query_hashes, both_strands=True)
    with open(output_file, "w") as fh:
        for i in xrange(0, len(query_sequence)-k + 1):
            kmer = query_sequence[i:i+k]
            p0 = 0
            p1 = 0
            p2 = 0
            rkmer = get_revcomp(kmer)
            if kmer in data:
                p1 = int(data[kmer])
            if rkmer in data:
                p2 = int(data[rkmer])
            p = max(p1,p2,p0)
            fh.write("%s\t%s\t%s\n" % (kmer, i, p))
   
def sc_compute_kmer_data(fasta_file, jellyfish_data_folder, jf_db, jf_dat, k, mintf, dumpmintf):
    '''
    '''
    ouput_prefix = os.path.join(
            jellyfish_data_folder,
            "%s___" % jf_db,
        )
    count_kmers(fasta_file, ouput_prefix, k, mintf=mintf)
    merge_kmers(jellyfish_data_folder, ouput_prefix, jf_db)
    dump_kmers(jf_db, jf_dat, dumpmintf=dumpmintf)
    print "Sort data..."
    temp_file = jf_dat+".temp"
    D = {
        "in": jf_dat,
        "out": temp_file,
    }
    command = "sort -k2nr %(in)s > %(out)s" % D
    print command
    os.system(command)
    command = "mv %(out)s %(in)s" % D
    print command
    os.system(command)
    