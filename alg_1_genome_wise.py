import numpy as np
from Bio import SeqIO 
import os
import time
import matplotlib.pyplot as plt

#Using Generated Datasets
sample_genomes = {}
for record in SeqIO.parse('Datasets\Spl_genomes_50_5\sample_genomes_1_50_s5.fasta', 'fasta'):
        sample_genomes[record.id] = str(record.seq)

#Getting biological datasets
def get_genomes(folder, file_format="fasta"):
    genome_dict = {}
    for file in os.listdir(folder):
        file_path = os.path.join(folder, file)
        if os.path.isfile(file_path):
            for record in SeqIO.parse(file_path, file_format):
                genome_dict[record.id] = str(record.seq)[:10000]
                break
    return genome_dict

directory = "D:\IIT Madras\Course material\8th Sem\J CS6024 Alg. App to Comp Bio\Project\Datasets\yeast_datasets"
genomes = get_genomes(directory, file_format="fasta")
print(len(genomes))


def get_neighborhood(kmer):
    nt_list = ['A','C','G','T']
    output = []
    for i in range(len(kmer)):
        nt_list = ['A','C','G','T']
        nt_list.remove(kmer[i])
        neighbors = [kmer[:i] + char + kmer[i+1:] for char in nt_list]
    
        output += neighbors
    return output

def extract_kmers(genome, k):
    unique_kmers = []
    for i in range(len(genome) - k + 1):
        kmer = genome[i:i+k]
        if 'N' in kmer:
            continue
        unique_kmers.append(kmer)
    return list(set(unique_kmers))

def get_rev_com_kmers(kmer_list):
    output = []
    nt_comp = {'A': 'T', 
               'C': 'G',
               'G': 'C',
               'T': 'A'}
    
    for kmer in kmer_list:
        rev_com = ''.join([nt_comp[i] for i in kmer])[::-1]
        output.append(rev_com)
    # print(len(kmer_list))
    output = list(set(output+kmer_list))
    # print(len(output))
    return output

def get_su_kmers(new_kmer_list, weak_kmer_list,su_kmer_dict):
    n = str(len(su_kmer_dict) + 1)
    su_kmers = []
    for kmer in new_kmer_list:
        flag = True
        if kmer in weak_kmer_list:
            continue
        else:
            flag = True
            for key in su_kmer_dict:
                if kmer in su_kmer_dict[key]:
                    su_kmer_dict[key].remove(kmer)
                    weak_kmer_list.append(kmer)
                    flag = False
                    break
        if not flag:
            continue
            
        flag = True
        neighbors = get_neighborhood(kmer)
        for neb in neighbors:
            if neb in weak_kmer_list:
                weak_kmer_list.append(kmer)
                flag = False
                
            for key in su_kmer_dict:
                if neb in su_kmer_dict[key]:
                    su_kmer_dict[key].remove(neb)
                    weak_kmer_list.append(neb)
                    flag = False
        if flag:
            su_kmers.append(kmer)
        
        
    su_kmer_dict['Genome_' + str(n)] = su_kmers

    return weak_kmer_list, su_kmer_dict

def run_full_alg_1(input, k):
    weak_kmers = []
    su_kmers = {}
    i = 1
    for key in input:
        print(i,key)
        genome_kmers = extract_kmers(input[key], k)
        full_kmer_list = get_rev_com_kmers(genome_kmers)
        weak_kmers, su_kmers = get_su_kmers(full_kmer_list, weak_kmers, su_kmers)
        i+=1
    for key in su_kmers:
        print(f'{key} : {len(su_kmers[key])}')

    return su_kmers

# start = time.time()
# strongly_unique_kmers = run_full_alg_1(genomes,31)
# end = time.time()
# print(end-start)


times = []
for i in [50,100,150,200,250]:
    sample_genomes = {}
    for record in SeqIO.parse(f'Datasets\diff_sim_datasets\sample_genomes_5k_{i}.fasta', 'fasta'):
            sample_genomes[record.id] = str(record.seq)
    
    start = time.time()
    strongly_unique_kmers = run_full_alg_1(sample_genomes,10)
    end = time.time()
    print(end-start)
    times.append(end-start)

print(times)

# Data for the plot





