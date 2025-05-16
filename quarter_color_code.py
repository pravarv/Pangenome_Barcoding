import numpy as np
from Bio import SeqIO 
import time
import os
import matplotlib.pyplot as plt
from collections import defaultdict
from collections import Counter

sample_genomes = {}
for record in SeqIO.parse('Datasets\Spl_genomes_50_5\sample_genomes_1_50_s5.fasta', 'fasta'):
        sample_genomes[record.id] = str(record.seq)

# genome_dict = {'Genome_1': 'ACAGATTC',
#                'Genome_2': 'ACAGATCC',
#                'Genome_3': 'AGAGATTC'}

def get_genomes(folder, file_format="fasta"):
    genome_dict = {}
    for file in os.listdir(folder):
        file_path = os.path.join(folder, file)
        if os.path.isfile(file_path):
            for record in SeqIO.parse(file_path, file_format):
                genome_dict[record.id] = str(record.seq)[:30000]
                break
    return genome_dict

directory = "D:\IIT Madras\Course material\8th Sem\J CS6024 Alg. App to Comp Bio\Project\Datasets\myco_plasm_genomes"
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

def pairwise_hamm_dist(kmer_list):
    strong_kmers = []
    weak_kmers = []
    k = len(kmer_list[0][0])
    for kmer in kmer_list:
        if all(sum([int(kmer[0][p] != other[0][p]) for p in range(k)]) >= 2 for other in kmer_list if kmer != other):
            strong_kmers.append(kmer)
        else:
            weak_kmers.append(kmer)

    return strong_kmers, weak_kmers

# print(pairwise_hamm_dist(['ACAGATTC','ACAGATCC','AGAGATTC']))

def extract_kmers(kmers, genome, k, n):

    for i in range(len(genome) - k + 1):
        kmer = genome[i:i+k]
        if 'N' in kmer:
            continue
        if kmer in kmers:
            kmers[kmer].append(n)
        else:
            kmers[kmer] = [n]
    return kmers

def get_rev_com_kmers(kmer_list):
    output = []
    nt_comp = {'A': 'T', 
               'C': 'G',
               'G': 'C',
               'T': 'A'}
    for kmer in kmer_list:
        rev_com = ''.join([nt_comp[i.upper()] for i in kmer[0]])[::-1]
        output.append((rev_com, kmer[1]))
    # print(len(kmer_list))
    output = output+kmer_list
    output = sorted(output)
    # print(len(output))
    return output

def rev_comp(kmer):
    nt_comp = {'A': 'T', 
               'C': 'G',
               'G': 'C',
               'T': 'A'}
    rev_com = ''.join([nt_comp[i.upper()] for i in kmer[0]])[::-1]
    return (rev_com,kmer[1])

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

def quarter_suk(kmer_list):
    # print('____')
    # print(len(kmer_list))
    kmer_list = kmer_list + get_rev_com_kmers(kmer_list)
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
        tribuckets[kmer[0][:(3*k//4)]].append(kmer)
    
    # print(len(tribuckets))
    for buck in tribuckets:
        if len(tribuckets[buck]) == 1:
            tri_su_kmers.append(tribuckets[buck][0])
        else:
            tri_suks,tri_weak = pairwise_hamm_dist(tribuckets[buck])
            tri_su_kmers += tri_suks
            tri_weak_kmers += tri_weak
    for kmer in tri_weak_kmers:
        if rev_comp(kmer)[0] in tri_su_kmers:
            tri_su_kmers.remove(rev_comp(kmer))
    
    # for key in tribuckets:
    #     print(f'{key} : {tribuckets[key]}')
    # print(len(tri_su_kmers))
    # print(tri_su_kmers)
    for kmer in kmer_list:
        quarbuckets[kmer[0][:k//2] +'|'+ kmer[0][(3*k//4):]].append(kmer)
    # print(len(quarbuckets))
    for buck in quarbuckets:
        if len(quarbuckets[buck]) == 1:
            quar_su_kmers.append(quarbuckets[buck][0])
        else:
            quar_suks, quar_weak = pairwise_hamm_dist(quarbuckets[buck])
            quar_su_kmers += quar_suks
            quar_weak_kmers += quar_weak
    for kmer in quar_weak_kmers:
        if rev_comp(kmer)[0] in quar_su_kmers:
            quar_su_kmers.remove(rev_comp(kmer)[0])
    # for key in quarbuckets:
    #     print(f'{key} : {quarbuckets[key]}')
    # print(len(quar_su_kmers))
    # print(quar_su_kmers)

    su_kmers = set(tri_su_kmers).intersection(set(quar_su_kmers))
    # print(len(su_kmers))
    # print(su_kmers)

    return su_kmers

# ukmers = []
# ukmers += extract_kmers(sample_genomes['Genome_6'], 5)
# su_kmers = quarter_suk(ukmers)
# print(sorted(su_kmers))

def run_full_alg_quarter(input, k):
    str_kmers = []
    su_kmers = {'Genome_' + str(p+1): [] for p in range(len(input))}
    
    pangenome_kmers = {}
    i = 1
    for key in input:
        print(i,key)
        pangenome_kmers = extract_kmers(pangenome_kmers,input[key], k, i)
        i += 1
    print(len(pangenome_kmers))
    weak_keys = []
    for key in pangenome_kmers:
        if len(pangenome_kmers[key]) > 1:
            weak_keys.append(key)
    
    for kmer in weak_keys:
        del pangenome_kmers[kmer]
    print(len(pangenome_kmers))
    
    all_kmers = [(kmer,pangenome_kmers[kmer][0]) for kmer in pangenome_kmers]
    # print(all_kmers[:5])

    str_kmers = quarter_suk(all_kmers)
    for kmer in str_kmers:
        su_kmers['Genome_' + str(kmer[1])].append(kmer[0])

    for key in su_kmers:
        print(f'{key} : {len(su_kmers[key])}')
    
    return su_kmers

start = time.time()
strongly_unique_kmers = run_full_alg_quarter(genomes,31)
end = time.time()
print(end-start)

# times = []
# for i in [50,100,150,200,250]:
#     sample_genomes = {}
#     for record in SeqIO.parse(f'Datasets\diff_sim_datasets\sample_genomes_5k_{i}.fasta', 'fasta'):
#             sample_genomes[record.id] = str(record.seq)
    
#     start = time.time()
#     strongly_unique_kmers = run_full_alg_quarter(sample_genomes,10)
#     end = time.time()
#     print(end-start)
#     times.append(end-start)

# print(times)







