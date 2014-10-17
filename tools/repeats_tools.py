#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 14.02.2014
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
"""
Functions related to repeats annotation in genomic data
"""

import os
import re
from trseeker.tools.repbase_tools import compute_repbase_stats
from trseeker.seqio.tab_file import sc_iter_tab_file
from trseeker.seqio.tab_file import sc_write_model_to_tab_file
from trseeker.models.annotation_models import BEDModel, GFF3Model


def run_windowmasker(input_fasta, output_dir):
    ''' Wrapper for running WindowMasker.

    windowmasker [-h] [-help] [-xmlhelp] [-ustat unit_counts]
    [-in input_file_name] [-out output_file_name] [-checkdup check_duplicates]
    [-fa_list input_is_a_list] [-mem available_memory] [-unit unit_length]
    [-genome_size genome_size] [-window window_size] [-t_extend T_extend]
    [-t_thres T_threshold] [-set_t_high score_value] [-set_t_low score_value]
    [-parse_seqids] [-outfmt output_format] [-t_high T_high] [-t_low T_low]
    [-infmt input_format] [-exclude_ids exclude_id_list] [-ids id_list]
    [-text_match text_match_ids] [-use_ba use_bit_array_optimization]
    [-sformat unit_counts_format] [-smem available_memory] [-dust use_dust]
    [-dust_level dust_level] [-mk_counts] [-convert] [-version-full]
    '''
    data = {
        "input_fasta": input_fasta,
        "output_mk_counts": os.path.join(output_dir, "windowmasker.mk_counts"),
        "output_wm": os.path.join(output_dir, "windowmasker.txt"),
        "output_bed": os.path.join(output_dir, "windowmasker.bed"),
        "output_gff3": os.path.join(output_dir, "windowmasker.gff3"),
    }
    command = "windowmasker -mk_counts -in %(input_fasta)s -out %(output_mk_counts)s" % data
    print command
    os.system(command)
    command = "windowmasker -in %(input_fasta)s -ustat %(output_mk_counts)s -outfmt interval -out %(output_wm)s" % data
    print command
    os.system(command)
    command = """awk 'BEGIN { OFS="\t"; scaffold="" } { if (match($0, ">.*")) { gsub(">", "", $0); gsub(" .*", "", $0); scaffold=$0 } else { gsub(" - ", "\t", $0); print scaffold "\t" $1 "\t" $2+1 } }' %(output_wm)s> %(output_bed)s""" % data
    print command
    os.system(command)

    gff_objs = []
    kwargs = {
        "source": "WindowMasker",
        "ftype": "Repeat",
        "score": 0,
        "strand": "+",
        "phase": "."

    }
    print "Convert to GFF3..."
    masked_bp = 0
    for bed_obj in sc_iter_tab_file(data["output_bed"], BEDModel):
        gff_obj = GFF3Model()
        gff_obj.set_with_bed_obj(bed_obj, **kwargs)
        gff_objs.append(gff_obj)
        masked_bp += abs(bed_obj.start - bed_obj.end)
    sc_write_model_to_tab_file(data["output_gff3"], gff_objs)
    return data["output_bed"], data["output_gff3"], masked_bp


def run_dustmasker(input_file, output_dir):
    ''' Wrapper for running dust.
    '''
    data = {
        'input_file': input_file,
        'output_file': os.path.join(output_dir, "dustmasker.txt"),
        'output_bed': os.path.join(output_dir, "dustmasker.bed"),
        'output_gff3': os.path.join(output_dir, "dustmasker.gff3"),
    }
    command = "dustmasker -in %(input_file)s -out %(output_file)s" % data
    print command
    os.system(command)
    command = """awk 'BEGIN { OFS="\t"; scaffold="" } { if (match($0, ">.*")) { gsub(">", "", $0); gsub(" .*", "", $0); scaffold=$0 } else { gsub(" - ", "\t", $0); print scaffold "\t" $1 "\t" $2+1 } }' %(output_file)s > %(output_bed)s""" % data
    print command
    print os.system(command)

    gff_objs = []
    kwargs = {
        "source": "DUST",
        "ftype": "SimpleSequence",
        "score": 0,
        "strand": "+",
        "phase": "."

    }
    print "Convert to GFF3..."
    masked_bp = 0
    for bed_obj in sc_iter_tab_file(data["output_bed"], BEDModel):
        gff_obj = GFF3Model()
        gff_obj.set_with_bed_obj(bed_obj, **kwargs)
        gff_objs.append(gff_obj)
        masked_bp += abs(bed_obj.start - bed_obj.end)
    sc_write_model_to_tab_file(data["output_gff3"], gff_objs)
    return data["output_bed"], data["output_gff3"], masked_bp


def run_repeatmasker_repbase(input_fastas, output_dir, species):
    ''' Wrapper for running RepeatMasker with Repbase database.
    TODO: remove hardcoded RM path.
    '''
    print "Running RepeatMasker for", species
    data = {
        "threads": 20,
        "species": species,
        "output_dir": output_dir,
        "input_fastas": input_fastas,
    }
    command = "/home/akomissarov/libs/RepeatMasker/RepeatMasker -s -pa %(threads)s -a -inv -xsmall -gff  -species %(species)s -dir %(output_dir)s  -nolow -source -html -u %(input_fastas)s" % data
    print command
    os.system(command)


def run_repeatmasker_local(input_fastas, output_dir, library_fasta):
    ''' Wrapper for running RepeatMasker with given database.
    TODO: remove hardcoded RM path.
    '''
    data = {
        "threads": 20,
        "library_fasta": library_fasta,
        "output_dir": output_dir,
        "input_fastas": input_fastas,
    }
    command = "/home/akomissarov/libs/RepeatMasker/RepeatMasker -s -pa %(threads)s -a -inv -xsmall -gff  -lib %(library_fasta)s -dir %(output_dir)s  -nolow -source -html -u %(input_fastas)s" % data
    print command
    os.system(command)









def run_ltr_harvest(input_fasta, output_dir):
    '''
    '''
    data = {
        "input_fasta": input_fasta,
        "index_name": os.path.join(output_dir, "indexname.fa"),
        "pred_out_file": os.path.join(output_dir, "ltrharvest.out.fa"),
        "pred_out_inner_file": os.path.join(output_dir, "ltrharvest.out.inner.fa"),
    }
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    if not os.path.isfile(pred_out_file):
        command = "/home/akomissarov/libs/genometools-unstable/bin/gt suffixerator -db %(input_fasta)s -indexname %(index_name)s -tis -suf -lcp -des -ssp -sds -dna" % data
        print command
        os.system(command)
        command = "/home/akomissarov/libs/genometools-unstable/bin/gt ltrharvest -index %(index_name)s -v -out %(pred_out_file)s -outinner %(pred_out_inner_file)s -seed 30 -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -minlenltr 100 -maxlenltr 1000 -mindistltr 1000 -maxdistltr 15000 -similar 90.0 -overlaps all -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -longoutput" % data
        print command
        os.system(command)
    command = "/home/akomissarov/libs/vmatch-2.2.3-Linux_x86_64-64bit/mkvtree -db %(pred_out_file)s -dna -pl -allout -v" % data
    print command
    os.system(command)

    command = "/home/akomissarov/libs/vmatch-2.2.3-Linux_x86_64-64bit/vmatch -dbcluster 95 7 Cluster-pred-chrAll -p -d -seedlength 50 -l 1101 -exdrop 9 %(pred_out_file)s" % data
    print command
    os.system(command)    






def run_repeatscout(input_fasta, output_dir):
    '''
    '''
    print input_fasta
    
    file_name = os.path.split(input_fasta)[-1]


    result = {}

    data = {
        "input_fasta": input_fasta,
        "output_freq": os.path.join(output_dir, "scout.freq"),
        "library_file": os.path.join(output_dir, "scout.lib.fa"),
        "filtered_library": os.path.join(output_dir, "scout.lib.f1.fa"),
        "repeat_masker_output_dir": os.path.join(output_dir, "repeat_masker"),
        "repeat_masker_not_trf_output_dir": os.path.join(output_dir, "repeat_masker_no_trf"),
        "repeatmasker_out_file": os.path.join(output_dir, "repeat_masker", "%s.tbl" % file_name),
        "rm_thresh": 10,
    }
    
    if not os.path.isfile(data["output_freq"]):
        command = "/home/akomissarov/libs/RepeatScout-1/build_lmer_table -sequence %(input_fasta)s -freq %(output_freq)s -v" % data
        print command
        os.system(command)

    if not os.path.isfile(data["library_file"]):
        command = "/home/akomissarov/libs/RepeatScout-1/RepeatScout -sequence %(input_fasta)s -output %(library_file)s -freq %(output_freq)s -vvvv" % data
        print command
        os.system(command)

    if not os.path.isfile(data["filtered_library"]):
        command = "cat %(library_file)s |  /home/akomissarov/libs/RepeatScout-1/filter-stage-1.prl > %(filtered_library)s" % data
        print command
        os.system(command)

    if not os.path.isdir(data["repeat_masker_output_dir"]):
        os.makedirs(data["repeat_masker_output_dir"])
        run_repeatmasker_local(data["input_fasta"], data["repeat_masker_output_dir"], data["filtered_library"])

    print data["repeatmasker_out_file"]
    if os.path.isfile(data["repeatmasker_out_file"]):
        with open(data["repeatmasker_out_file"]) as fh:
            data = fh.read()
            masked_percent = re.findall("bases masked:\s*(\d+)\s*bp\s*\(\s*(\d+.\d+)\s*\%\)", data, re.S)
            print masked_percent
            result["masked"] = int(masked_percent[0][0])
            result["pmasked"] = float(masked_percent[0][1])
    else:
        result["masked"] = 0
        result["pmasked"] = 0.0





    # command = "cat %(filtered_library)s |  /home/akomissarov/libs/RepeatScout-1/filter-stage-2.prl --cat %(repeatmasker_out_file)s --thresh=%(rm_thresh)s" % data
    # print command
    # os.system(command)  

    # compile families
    

    # compare with known tandem repeats by blast

    # compare with known repeats by blast
    # for family_fasta_file in iter_folder()

    return result


def run_process_repeatscout_results(input_fasta, output_dir):
    '''
    '''
    only_file_name = os.path.split(input_fasta)[-1]
    data = {
        "repeatmasker_ori_out_file": os.path.join(output_dir, "repeat_masker", "%s.ori.out" % only_file_name),
        "output_families_folder": os.path.join(output_dir, "families"),
        "output_tab_file": os.path.join(output_dir, "scout.families.tsv"),
    }
    if not os.path.isdir(data["output_families_folder"]):
        os.makedirs(data["output_families_folder"])

    
    compute_repbase_stats(data["repeatmasker_ori_out_file"], data["output_tab_file"], input_fasta, data["output_families_folder"])
    

def run_repeat_modeler(db_name, input_file):
    '''
    '''
    data = {
        'db_name': db_name,
        'input_file': input_file,
    }
    command = "/home/akomissarov/libs/RepeatModeler/BuildDatabase -name %(db_name)s -engine ncbi %(input_file)s" % data
    print command
    os.system(command)
    command = "/home/akomissarov/libs/RepeatModeler/RepeatModeler -database %(db_name)s -engine ncbi -pa 10" % data
    print command
    print os.system(command)






