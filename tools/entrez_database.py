#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 11.02.2022
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com
'''
### Function for downloading datasets from NCBI

```python
from trseeker.tools.entrez_datanase import *
```

Download all proteins according to a given query into output file and then return these proteins as fasta text.

```python
fasta_data = download_proteins_from_ncbi(query, output_file, email, batch=500, verbose_step=1000)
```

Download items from given database according to query, rettype, and retmode. Save output to output_file and return as text.

```python
data = download_items_from_ncbi(query, 
                             database, 
                             output_file, 
                             email, 
                             rettype="fasta", 
                             retmode="text", 
                             batch=500, 
                             verbose_step=1000)

```

Get items from given database according to query, rettype, and retmode.

```python
items = get_items_from_ncbi(query, 
                         database, 
                         email, 
                         rettype="fasta", 
                         retmode="text", 
                         batch=500, 
                         verbose_step=1000)
```

Get RNA SRA datasets from NCBI according to taxid. 

Output includes LIBRARY_SOURCE, STUDY_ABSTRACT, DESIGN_DESCRIPTION, PRIMARY_ID, DESCRIPTION, LINKS.

And the second part of the output contains full xml data.

```python
data, _ = get_rna_sra_datasets_by_taxid(taxid, email, batch=500, verbose_step=1)

all_links = {}

for *item, links in datasets.values():
    links = dict([
                   (url.split("/")[-1], url)
                        for 
                            url 
                                in links])
    for sra in links:
        all_links[sra] = links[sra]
```

Download and unpack SRA filrs from NCBI according to taxid.

```python
download_rna_sra_datasets_by_taxid(taxid, email, output_folder, threads=30, batch=500, verbose_step=1)
```

Download genomes and annotation from NCBI according to taxid.
Return refseq and genbank datasets.

```python
download_genome_assemblies_and_annotation_from_ncbi(taxid, output_folder, threads=30, only_refseq=True)
```

'''

from Bio import Entrez
import os
import tempfile
import shutil
import urllib.request
import re
import time
from urllib.error import HTTPError


class NoToolException(Exception):
    pass


def _get_feature(feature, record):
    d = re.findall("<%s>(.*?)</%s>" % (feature, feature), record, re.M)
    if d:
        return d[0]
    return ""


def _download_genomic_links(url, only_gff=False, add_fasta=True):
    '''
    '''
    if not url.endswith("/"):
        url += "/"
    url = url.replace("ftp://", " https://")
    to_download = []
    print(f"Parsing {url} ...", end=" ")
    with urllib.request.urlopen(url) as response:
        html = response.read().decode("utf8")
        links = re.findall("<a href=\"(.*?)\">(.*?)</a>", html, re.M|re.S)
        for link in links:
            if only_gff:
                if "_genomic.gff.gz" in link[0]:
                    to_download.append(url + link[0])    
                elif add_fasta and  "_genomic.fna.gz" in link[0]:
                    to_download.append(url + link[0])
                else:
                    continue
            if "_genomic.gff.gz" in link[0]:
                to_download.append(url + link[0])
            elif "_genomic.fna.gz" in link[0]:
                to_download.append(url + link[0])
            elif "_protein.faa.gz" in link[0]:
                to_download.append(url + link[0])
            elif "translated_cds.faa.gz" in link[0]:
                to_download.append(url + link[0])
    print(f" found {len(to_download)} links.")
    return to_download


def is_tool(program_name):
    '''Check whether program_name is on PATH and executable.
    '''
    return shutil.which(program_name) is not None


def download_proteins_from_ncbi(query, output_file, email, batch=500, verbose_step=1000):
    ''' Download all proteins according to a given query into output file 
        and then return these proteins as fasta text.
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
                             rettype="fasta", 
                             retmode="text", 
                             batch=500, 
                             verbose_step=1000):
    ''' Download items from given database according to query, rettype, and retmode. 
        Save output to output_file and return as text.
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
                         rettype="fasta", 
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
    ''' Get RNA SRA datasets from NCBI. 
    Output includes 
        LIBRARY_SOURCE, 
        STUDY_ABSTRACT, 
        DESIGN_DESCRIPTION, 
        PRIMARY_ID, 
        DESCRIPTION, 
        LINKS
    And the second part of the output contains full xml data.
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


def download_rna_sra_datasets_by_taxid(taxid, email, output_folder, threads=30, batch=500, verbose_step=1, quiet=True):
    ''' Download and unpack SRA filrs from NCBI according to taxid.
    '''
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)

    if not is_tool("fastq-dump"):
        raise NoToolException("Please, install 'mamba install -c bioconda sra-tools'")

    data, _ = get_rna_sra_datasets_by_taxid(taxid, email, batch=500, verbose_step=1)
    all_links = []
    for *item, links in data.values():
        links = dict([
                       (url.split("/")[-1], url)
                            for 
                                url 
                                    in links])
        all_links.append(links)

    file_with_link = os.path.join(output_folder, "to_download.list")
    with open(file_with_link, "w") as fw:
        for links in all_links:
            for dataset, url in links.items():
                fw.write(f"{url}\n")

    os.chdir(output_folder)
    if quiet:
        command = f"less to_download.list | xargs -P {threads} -n 1 wget -q"
    else:
        command = f"less to_download.list | xargs -P {threads} -n 1 wget"
    print(command)
    os.system(command)

    os.chdir(output_folder)
    command = f"less to_download.list | rev | cut -f1 -d"/" | rev | xargs -P {threads} -n 1 fastq-dump --split-files"
    print(command)
    os.system(command)


def _is_file_exists(url, output_folder):
    ''' Check file name existence.'''
    file_name = url.split("/")[-1]
    output_file = os.path.join(output_folder, file_name)
    output_file_unpacked = output_file
    if output_file.endswith(".gz"):
        output_file_unpacked = output_file.replace(".gz", "")
    if os.path.isfile(output_file) or os.path.isfile(output_file_unpacked):
        print(f"File exists: {file_name}")
        return True
    return False


def download_genome_assemblies_and_annotation_from_ncbi(taxid, output_folder, threads=30, only_refseq=True, only_gff=False, add_fasta=False, quiet=True, mock=False, assemblies_to_use=None):
    ''' Download genomes and annotation from NCBI according to taxid.
        Return refseq and genbank datasets.
    '''
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
        
    refseq_results = {}
    genbank_results = {}
    errors = []
    with tempfile.NamedTemporaryFile() as temp_output_file:
        if not is_tool("esearch"):
            raise NoToolException("Please, install 'mamba install -c bioconda entrez-direct'")
        
        command = f"esearch -db assembly -query '{taxid}[Organism:exp]' | efetch -format docsum > {temp_output_file.name}"
        print(command)
        os.system(command)
        
        temp_output_file.seek(0)
        data = temp_output_file.read().decode("utf8")

        documents = re.findall("<DocumentSummary>(.*?)</DocumentSummary>", data, re.S)

        print(f"Found {len(documents)} items...")
        
        for i, document in enumerate(documents):
            print(f"Progress {i}/{len(documents)} found {len(refseq_results)} RefSeq links and {len(genbank_results)} GenBank links and {len(errors)} errors.")

            
            genbank_paths = re.findall('<Taxid>(\d+)</Taxid>.*?<Organism>(.*?)</Organism>.*?<FtpPath_GenBank>(.*?)</FtpPath_GenBank>', document, re.S)
            refseq_paths = re.findall('<Taxid>(\d+)</Taxid>.*?<Organism>(.*?)</Organism>.*?<FtpPath_RefSeq>(.*?)</FtpPath_RefSeq>', document, re.S)
                        
            for taxid, organism, url in refseq_paths:
                assembly = url.split("/")[-1]
                if assemblies_to_use and not assembly in assemblies_to_use:
                    continue
                name = (taxid, organism, assembly)
                while True:
                    try:
                        refseq_results[name] = _download_genomic_links(url, only_gff=only_gff, add_fasta=add_fasta)
                        break
                    except HTTPError as err:
                        if err.code == 404:
                           errors.append(url)
                           break
                        else:
                           raise err
                    except:
                        print("Retry...")
                        time.sleep(5)
                
                
            if not only_refseq:
                for taxid, organism, url in genbank_paths:
                    assembly = url.split("/")[-1]
                    if assemblies_to_use and not assembly in assemblies_to_use:
                        continue
                    name = (taxid, organism, assembly)
                    while True:
                        try:
                            genbank_results[name] = _download_genomic_links(url, only_gff=only_gff, add_fasta=add_fasta)
                            break
                        except HTTPError as err:
                            if err.code == 404:
                               errors.append(url)
                               break
                            else:
                               raise err
                        except:
                            print("Retry...")
                            time.sleep(5)

        
        print(f"Found {len(refseq_results)} RefSeq links and {len(genbank_results)} GenBank links.")

    file_with_link = os.path.join(output_folder, "to_download.list")
    links_to_download = 0
    with open(file_with_link, "w") as fw:
        for organism in refseq_results:
            for url in refseq_results[organism]:
                if _is_file_exists(url, output_folder):
                    continue
                fw.write(f"{url}\n")
                links_to_download += 1
                
        if not only_refseq:
            for organism in genbank_results:
                if organism in refseq_results:
                    continue
                for url in genbank_results[organism]:
                    if _is_file_exists(url, output_folder):
                        continue
                    fw.write(f"{url}\n")
                    links_to_download += 1

    if mock or links_to_download == 0:
        return refseq_results, genbank_results
    print(f"Downloading {links_to_download} links")
    os.chdir(output_folder)
    if quiet:
        command = f"less to_download.list | xargs -P {threads} -n 1 wget -q"
    else:
        command = f"less to_download.list | xargs -P {threads} -n 1 wget"
    print(command)
    os.system(command)
    
    print("errors:")
    for error in errors:
        print(error)

    return refseq_results, genbank_results