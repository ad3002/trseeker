'''
Jellyfish python wrapper
'''
import os
from trseeker.settings import load_settings
from PyExp.readers.abstract_reader import sc_iter_filepath_folder
import subprocess

settings = load_settings()
location = settings["blast_settings"]["jellyfish_location"]

def count_kmers(input_fasta, ouput_prefix, k, mintf=None):
    ''' Count kmers with Jellyfish
    '''
    params = {
        "location": location,
        "input_fasta": input_fasta,
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
    params = {
        "location": location,
        "ouput_file": ouput_file+".jf",
        "ouput_prefix": ouput_prefix,
    }
    file_count = 0
    for file_name in sc_iter_filepath_folder(folder, mask="."):
        head, tail = os.path.split(ouput_file)
        if tail in file_name:
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
    params = {
        "location": location,
        "db_file": db_file,
        "stats_file": stats_file,
    }
    command = "%(location)s stats --verbose -o %(stats_file)s %(db_file)s" % params
    print "Execute:", command
    os.system(command)

def histo_kmers(db_file, histo_file):
    params = {
        "location": location,
        "db_file": db_file,
        "histo_file": histo_file,
    }
    command = "%(location)s histo --verbose -o %(histo_file)s %(db_file)s" % params
    print "Execute:", command
    os.system(command)

def dump_kmers(db_file, fasta_file):
    params = {
        "location": location,
        "db_file": db_file,
        "fasta_file": fasta_file,
    }
    command = "%(location)s dump --column --tab -o %(fasta_file)s %(db_file)s" % params
    print "Execute:", command
    os.system(command)

def query_kmers(db_file, query_hashes, both_strands=True):
    '''
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