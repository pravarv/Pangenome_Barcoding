import numpy as np
from Bio import SeqIO 
import os
import time
import matplotlib.pyplot as plt

#Using Generated Datasets
sample_genomes = {}
for record in SeqIO.parse('Datasets\Spl_genomes_2k_20\sample_genomes_1_2k_s20.fasta', 'fasta'):
        sample_genomes[record.id] = str(record.seq)

#Getting biological datasets
def get_genomes(folder, file_format="fasta"):
    genome_dict = {}
    for file in os.listdir(folder):
        file_path = os.path.join(folder, file)
        if os.path.isfile(file_path):
            for record in SeqIO.parse(file_path, file_format):
                genome_dict[record.id] = str(record.seq)[:20000]
                break
    return genome_dict

directory = "D:\IIT Madras\Course material\8th Sem\J CS6024 Alg. App to Comp Bio\Project\Datasets\Yeast_genomes"
genomes = get_genomes(directory, file_format="fasta")
# for key in genomes:
#     print(f'{key}: {len(genomes[key])}')

#Required Functions

def get_neighborhood(kmer):
    nt_list = ['A','C','G','T']
    output = []
    for i in range(len(kmer)):
        nt_list = ['A','C','G','T']
        nt_list.remove(kmer[i].upper())
        neighbors = [kmer[:i] + char + kmer[i+1:] for char in nt_list]
        output += neighbors
    return output

def rev_comp(kmer):
    nt_comp = {'A': 'T', 
               'C': 'G',
               'G': 'C',
               'T': 'A'}
    rev_com = ''.join([nt_comp[i.upper()] for i in kmer])[::-1]
    return rev_com

def extract_kmers(genome, k):
    unique_kmers = []
    for i in range(len(genome) - k + 1):
        kmer = genome[i:i+k]
        if 'N' in kmer:
            continue
        unique_kmers.append(sorted([kmer,rev_comp(kmer)])[0])
    unique_kmers = list(set(unique_kmers))
    # print('Length: ', len(unique_kmers))
    return unique_kmers

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

def get_kmer_index(input_genomes, k):
    genome_kmers = extract_kmers(input_genomes[sorted(list(input_genomes.keys()))[0]], k)
    # full_kmer_list = get_rev_com_kmers(genome_kmers)
    kmers = np.array(genome_kmers)
    kmers = np.sort(kmers)
    # kmers = np.array(full_kmer_list)
    index = np.ones((len(kmers),1))
    for key in input_genomes:
        print(key)
        if key == sorted(list(input_genomes.keys()))[0]:
            continue
        genome_kmers = extract_kmers(input_genomes[key], k)
        # full_kmer_list = get_rev_com_kmers(genome_kmers)
        index = np.hstack((index,np.zeros((kmers.shape[0],1))))
        i = 1
        for kmer in genome_kmers:
            b = np.searchsorted(kmers, kmer)
            if b == kmers.shape[0]:
                kmers = np.concatenate((kmers, np.array([kmer])))
                index = np.vstack((index, np.array([0]*(index.shape[1]-1)+[1])))
            else:
                if kmers[b] == kmer:
                    index[b,index.shape[1]-1] = 1
                else:
                    kmers = np.concatenate((kmers[:b], np.array([kmer]),kmers[b:]))
                    index = np.vstack((index[:b], np.array([0]*(index.shape[1]-1)+[1]), index[b:]))
            i+=1
            if i % 10000 == 0:
                print(i)

    # print(kmers,index)
    return kmers, index

def get_unique_kmers(kmers, index):
    unique_kmer_dict = {'Genome_' + str(i+1) : [] for i in range(index.shape[1])}
    for i in range(index.shape[0]):
        if np.sum(index[i]) == 1.0:
            x = np.where(index[i] == 1.0)
            unique_kmer_dict['Genome_' + str(int(x[0]+1))].append(kmers[i])
                
    return unique_kmer_dict

def get_strongly_unique_kmers(kmer_dict, kmers):
    for key in kmer_dict:
        flagged = []
        for kmer in kmer_dict[key]:
            nbrs = get_neighborhood(kmer)
            for neb in nbrs:
                if neb in kmers:
                    flagged.append(kmer)
                if rev_comp(neb) in kmers:
                    flagged.append(kmer)
        for x in set(flagged):
            kmer_dict[key].remove(x)
    
    return kmer_dict

#Main Function: Input: Genomes,k ; Output: Strongly Unique Kmers in each Genome
def run_full_alg_2(input,k):
    # for key in input:
    #     print(f'{key} : {input[key]}')
    kmer_array, kmer_index = get_kmer_index(input,k)
    unique_kmers = get_unique_kmers(kmer_array,kmer_index)
    # for key in unique_kmers:
    #     print(f'{key} : {len(unique_kmers[key])})

    su_kmers = get_strongly_unique_kmers(unique_kmers,kmer_array)
    for key in su_kmers:
        print(f'{key},  : {len(su_kmers[key])}')

    return su_kmers

# start = time.time()
# strongly_unique_kmers = run_full_alg_2(genomes,31)
# end = time.time()
# print(end-start)

# Generating Variable Similarity Runtimes
times = []
for i in [50,100,150,200,250]:
    sample_genomes = {}
    for record in SeqIO.parse(f'Datasets\diff_sim_datasets\sample_genomes_2k_{i}.fasta', 'fasta'):
            sample_genomes[record.id] = str(record.seq)
    
    start = time.time()
    strongly_unique_kmers = run_full_alg_2(sample_genomes,15)
    end = time.time()
    print(end-start)
    times.append(end-start)
print(times)

# Generating Variable Genome Size Runtimes
times = []
for i in [1,2,5,10,15]:
    sample_genomes = {}
    for record in SeqIO.parse(f'Datasets\diff_size_datasets\sample_genomes_{i}k_1.fasta', 'fasta'):
            sample_genomes[record.id] = str(record.seq)
    
    start = time.time()
    strongly_unique_kmers = run_full_alg_2(sample_genomes,15)
    end = time.time()
    print(end-start)
    times.append(end-start)
print(times)