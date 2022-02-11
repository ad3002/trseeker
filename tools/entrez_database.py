#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 11.02.2022
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com
'''
Downloading datasets from NCBI


- download_proteins_from_ncbi(query, output_file, email, batch=500, verbose_step=1000)


# https://www.ncbi.nlm.nih.gov/assembly/?term=txid3699[Organism:exp]
# https://linsalrob.github.io/ComputationalGenomicsManual/Databases/NCBI_Edirect.html

esearch -db assembly -query "txid3699[Organism:exp]" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F"/" '{print "curl -o "$NF"_genomic.fna.gz " $0"/"$NF"_genomic.fna.gz"}'
'''

from Bio import Entrez
import re
import tempfile

def _get_feature(feature, record):
    d = re.findall("<%s>(.*?)</%s>" % (feature, feature), record, re.M)
    if d:
        return d[0]
    return ""


def download_proteins_from_ncbi(query, output_file, email, batch=500, verbose_step=1000):
    '''
    '''
    Entrez.email = email
    
    start = 0
    proteins = []
    while True:
        handle = Entrez.esearch(db='protein', retmax=batch, retstart=start, term=query)
        record = Entrez.read(handle)
        handle.close()
        start += batch
        if start % verbose_step == 0:
            print(len(proteins), end=" ")
        proteins += record["IdList"]
        if len(record["IdList"]) == 0:
            break
    print()
    print(f"Downloaded {len(proteins)} proteins")

    with open(output_file, "w") as fw:
        for i in range(0,len(proteins), batch):
            if i % verbose_step == 0:
                print(i, end=" ")
            handle = Entrez.efetch(db="protein", 
                                   id=proteins[i:i+batch], 
                                   rettype="fasta", 
                                   retmode="text")
            record = handle.read().replace("\n\n","\n")
            fw.write(record.strip())
            fw.write("\n")
    return proteins


def download_items_from_ncbi(query, 
                             database, 
                             output_file, 
                             email, 
                             rettype="tasta", 
                             retmode="text", 
                             batch=500, 
                             verbose_step=1000):
    '''
    '''
    Entrez.email = email
    
    start = 0
    items = []
    while True:
        handle = Entrez.esearch(db=database, retmax=batch, retstart=start, term=query)
        record = Entrez.read(handle)
        handle.close()
        start += batch
        if start % verbose_step == 0:
            print(len(items), end=" ")
        items += record["IdList"]
        if len(record["IdList"]) == 0:
            break
    print()
    print(f"Downloaded {len(items)} items")

    with open(output_file, "w") as fw:
        for i in range(0,len(items), batch):
            if i % verbose_step == 0:
                print(i, end=" ")
            handle = Entrez.efetch(db=database, 
                                   id=items[i:i+batch], 
                                   rettype=rettype, 
                                   retmode=retmode)
            record = handle.read().replace("\n\n","\n")
            fw.write(record.strip())
            fw.write("\n")
    return items


def get_items_from_ncbi(query, 
                         database, 
                         email, 
                         rettype="tasta", 
                         retmode="text", 
                         batch=500, 
                         verbose_step=1000):
    '''
    '''
    Entrez.email = email
    
    start = 0
    items = []
    while True:
        handle = Entrez.esearch(db=database, retmax=batch, retstart=start, term=query)
        record = Entrez.read(handle)
        handle.close()
        start += batch
        if start % verbose_step == 0:
            print(len(items), end=" ")
        items += record["IdList"]
        if len(record["IdList"]) == 0:
            break
    print()
    print(f"Downloaded {len(items)} items")

    resutls = []
    for i in range(0,len(items), batch):
        if i % verbose_step == 0:
            print(i, end=" ")
        handle = Entrez.efetch(db=database, 
                               id=items[i:i+batch], 
                               rettype=rettype, 
                               retmode=retmode)
        record = handle.read().replace("\n\n","\n")
        resutls.append(record)
    return items, resutls



def get_rna_sra_datasets_by_taxid(taxid, 
                         email, 
                         batch=500, 
                         verbose_step=1):
    '''
    '''
    Entrez.email = email
    
    query = '%s[Organism:exp] AND ("biomol rna"[Properties])' % taxid

    start = 0
    datasets = []
    while True:
        handle = Entrez.esearch(db="sra", retmax=batch, retstart=start, term=query)
        record = Entrez.read(handle)
        handle.close()
        start += batch
        if start % verbose_step == 0:
            print(len(datasets), end=" ")
        datasets += record["IdList"]
        if len(record["IdList"]) == 0:
            break
    print()
    print(f"Downloaded {len(datasets)} datasets")

    ALL_RNA_DATASETS = {}
    ALL_RNA_DATASETS_FULL = {}

    for i, seq_id in enumerate(datasets):
        if i % verbose_step == 0:
            print(i, end=" ")
        handle = Entrez.efetch(db="sra", 
                               id=seq_id, 
                               rettype="None", 
                               retmode="xml")
        record = handle.read().decode("utf8")
        LIBRARY_SOURCE = _get_feature("LIBRARY_SOURCE", record)
        STUDY_ABSTRACT = _get_feature("STUDY_ABSTRACT", record)
        DESIGN_DESCRIPTION = _get_feature("DESIGN_DESCRIPTION", record)
        PRIMARY_ID = _get_feature("PRIMARY_ID", record)
        DESCRIPTION = _get_feature("DESCRIPTION", record)
        LINKS = set(re.findall('url=\"(https://sra-downloadb.*?)\"', record, re.M))
        if LIBRARY_SOURCE == 'GENOMIC':
            continue
        
        ALL_RNA_DATASETS_FULL[seq_id] = record
        ALL_RNA_DATASETS[seq_id] = [LIBRARY_SOURCE, STUDY_ABSTRACT, DESIGN_DESCRIPTION, PRIMARY_ID, DESCRIPTION, LINKS]

    return ALL_RNA_DATASETS, ALL_RNA_DATASETS_FULL


def download_genome_assemblies_and_annotation_from_ncbi(taxid, 
                         email, 
                         batch=500, 
                         verbose_step=1):
    '''
    '''
    temp_output_file = tempfile.NamedTemporaryFile(
    command = f"esearch -db assembly -query '{taxid}[Organism:exp]' | efetch -format docsum > {temp_output_file}"
