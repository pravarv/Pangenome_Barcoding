import numpy as np
from Bio import SeqIO 
import os
import time
import matplotlib.pyplot as plt
from collections import defaultdict

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
                genome_dict[record.id] = str(record.seq)[80000:100000]

                break
    return genome_dict
folder = "D:\IIT Madras\Course material\8th Sem\J CS6024 Alg. App to Comp Bio\Project\Datasets\Myco_genomes"
genomes = get_genomes(folder, file_format="fasta")
# for key in genomes:
#     print(f'{key}: {len(genomes[key])}')

#Required Functions

def pairwise_hamm_dist(kmer_list):
    strong_kmers = []
    weak_kmers = []
    k = len(kmer_list[0])
    for kmer in kmer_list:
        if all(sum([int(kmer[p] != other[p]) for p in range(k)]) >= 2 for other in kmer_list if kmer != other):
            strong_kmers.append(kmer)
        else:
            weak_kmers.append(kmer)

    return strong_kmers, weak_kmers

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
        rev_com = ''.join([nt_comp[i.upper()] for i in kmer])[::-1]
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
        # print(key)
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

def quarter_suk(kmer_list):
    # print('____')
    # print(len(kmer_list))
    kmer_list = list(set(kmer_list + get_rev_com_kmers(kmer_list)))
    kmer_list = sorted(kmer_list)
    # print(len(kmer_list))
    tribuckets = defaultdict(list)
    quarbuckets = defaultdict(list)
    tri_su_kmers = []
    tri_weak_kmers = []
    quar_su_kmers = []
    quar_weak_kmers = []
    n = len(kmer_list)
    k = len(kmer_list[0])
    for kmer in kmer_list:
        tribuckets[kmer[:(3*k//4)]].append(kmer)
    
    # print(len(tribuckets))
    for buck in tribuckets:
        if len(tribuckets[buck]) == 1:
            tri_su_kmers.append(tribuckets[buck][0])
        else:
            tri_suks,tri_weak = pairwise_hamm_dist(tribuckets[buck])
            tri_su_kmers += tri_suks
            tri_weak_kmers += tri_weak
    for kmer in tri_weak_kmers:
        if rev_comp(kmer) in tri_su_kmers:
            tri_su_kmers.remove(rev_comp(kmer))
    
    # for key in tribuckets:
    #     print(f'{key} : {tribuckets[key]}')
    # print(len(tri_su_kmers))
    # print(tri_su_kmers)
    for kmer in kmer_list:
        quarbuckets[kmer[:k//2] +'|'+ kmer[(3*k//4):]].append(kmer)
    # print(len(quarbuckets))
    for buck in quarbuckets:
        if len(quarbuckets[buck]) == 1:
            quar_su_kmers.append(quarbuckets[buck][0])
        else:
            quar_suks, quar_weak = pairwise_hamm_dist(quarbuckets[buck])
            quar_su_kmers += quar_suks
            quar_weak_kmers += quar_weak
    for kmer in quar_weak_kmers:
        if rev_comp(kmer) in quar_su_kmers:
            quar_su_kmers.remove(rev_comp(kmer))
    # for key in quarbuckets:
    #     print(f'{key} : {quarbuckets[key]}')
    # print(len(quar_su_kmers))
    # print(quar_su_kmers)

    su_kmers = set(tri_su_kmers).intersection(set(quar_su_kmers))
    # print(len(su_kmers))
    # print(su_kmers)

    return su_kmers

def get_strongly_unique_kmers(kmer_dict, full_suks):
    for key in kmer_dict:
        flagged = []
        for kmer in kmer_dict[key]:
            if kmer in full_suks:
                continue
            else:
                flagged.append(kmer)
        
        for x in set(flagged):
            kmer_dict[key].remove(x)
     
    return kmer_dict

#Main Function: Input: Genomes,k ; Output: Strongly Unique Kmers in each Genome
def run_full_alg_3(input,k):
    # for key in input:
    #     print(f'{key} : {input[key]}')
    kmer_array, kmer_index = get_kmer_index(input,k)
    unique_kmers = get_unique_kmers(kmer_array,kmer_index)
    all_suks = quarter_suk(kmer_array.tolist())
    # for key in unique_kmers:
    #     print(f'{key} : {len(unique_kmers[key])})

    su_kmers = get_strongly_unique_kmers(unique_kmers,all_suks)

    for key in su_kmers:
        print(f'{key},  : {len(su_kmers[key])}')

    return su_kmers

start = time.time()
strongly_unique_kmers = run_full_alg_3(genomes,31)
end = time.time()
print(end-start)

# Generating Variable Similarity Runtimes
# times = []
# for i in [50,100,150,200,250]:
#     sample_genomes = {}
#     for record in SeqIO.parse(f'Datasets\diff_sim_datasets\sample_genomes_2k_{i}.fasta', 'fasta'):
#             sample_genomes[record.id] = str(record.seq)
    
#     start = time.time()
#     strongly_unique_kmers = run_full_alg_3(sample_genomes,15)
#     end = time.time()
#     print(end-start)
#     times.append(end-start)
# print(times)

# Generating Variable Genome Size Runtimes
# times = []
# for i in [1,2,5,10,15]:
#     sample_genomes = {}
#     for record in SeqIO.parse(f'Datasets\diff_size_datasets\sample_genomes_{i}k_1.fasta', 'fasta'):
#             sample_genomes[record.id] = str(record.seq)
    
#     start = time.time()
#     strongly_unique_kmers = run_full_alg_3(sample_genomes,15)
#     end = time.time()
#     print(end-start)
#     times.append(end-start)
# print(times)

# k_dict = {}
# for lenk in range(15,25):
#     start = time.time()
#     strongly_unique_kmers = run_full_alg_3(genomes,lenk)
#     k_dict[lenk] = [len(strongly_unique_kmers[key]) for key in strongly_unique_kmers]
#     end = time.time()
#     print(end-start)

# plot_list = defaultdict(list)
# for key in k_dict:
#     for i in range(5):
#         plot_list[i+1].append(k_dict[key][i])

# for key in plot_list:
#     print(f'{key}: {plot_list[key]}')








