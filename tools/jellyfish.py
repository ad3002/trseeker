'''
Jellyfish python wrapper
'''
import os
from trseeker.settings import load_settings
from PyExp import sc_iter_filepath_folder
import subprocess
from trseeker.tools.seqfile import sort_file_by_int_field

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
        "hash_size": 30000000,
        "hash_bits": 3,
        "threads": 8,
        "both_strands": "--both-strands",
        "ouput_prefix": ouput_prefix,
        "mintf": "",
    }
    if mintf:
        params["mintf"] = "--lower-count=%s" % mintf
    command = "%(location)s count %(mintf)s -m %(k)s -o %(ouput_prefix)s -c %(hash_bits)s -s %(hash_size)s %(both_strands)s -t %(threads)s %(input_fasta)s" % params
    print "Execute:", command
    os.system(command)

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

    # command = "rm %(ouput_prefix)s_*" % params
    # os.system(command)

    # return params["ouput_file"]

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

def dump_kmers(db_file, kmers_file):
    '''Dump k-mers database to tab-delimited file.
    '''
    params = {
        "location": location,
        "db_file": db_file,
        "kmers_file": kmers_file,
    }
    command = "%(location)s dump --column --tab -o %(kmers_file)s %(db_file)s" % params
    print "Execute:", command
    os.system(command)
    sort_file_by_int_field(kmers_file, 1)

def query_kmers(db_file, query_hashes, both_strands=True):
    '''Query jellyfish database.
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
                key, value = item.strip().split()
                final_result[key] = value
            else:
                final_result[key] = data[1]
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
        p = 0
        for i in xrange(0, len(DD)-18):
            kmer = query_sequence[i:i+17]
            if kmer in data:
                if int(data[kmer]) > 0:
                    p += 2
            p -= 1
            if p < 0:
                p = 0
            fh.write("%s\t%s\n" % (i, "*"*p))
   