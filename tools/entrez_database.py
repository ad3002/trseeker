#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 11.02.2022
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com
'''
Downloading datasets from NCBI


- download_proteins_from_ncbi(query, output_file, email, batch=500, verbose_step=1000)

'''

from Bio import Entrez


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
                print(len(i), end=" ")
            handle = Entrez.efetch(db="protein", 
                                   id=proteins[i:i+batch], 
                                   rettype="fasta", 
                                   retmode="text")
            record = handle.read().replace("\n\n","\n")
            fw.write(record.strip())
            fw.write("\n")
    return proteins




