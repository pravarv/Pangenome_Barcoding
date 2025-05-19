import numpy as np
from Bio import SeqIO 
import time
import os
import matplotlib.pyplot as plt
from collections import defaultdict
from collections import Counter

sample_genomes = {}
for record in SeqIO.parse('Datasets/Spl_genomes_50_5/test.fasta', 'fasta'):
        sample_genomes[record.id] = str(record.seq)

# def get_genomes(folder, file_format="fasta"):
#     genome_dict = {}
#     for file in os.listdir(folder):
#         file_path = os.path.join(folder, file)
#         if os.path.isfile(file_path):
#             for record in SeqIO.parse(file_path, file_format):
#                 genome_dict[record.id] = str(record.seq)[:10000]
#                 break
#     return genome_dict

# directory = "D:\IIT Madras\Course material\8th Sem\J CS6024 Alg. App to Comp Bio\Project\Datasets\yeast_datasets"
# genomes = get_genomes(directory, file_format="fasta")
# print(len(genomes))

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
    k = len(kmer_list[0])
    for kmer in kmer_list:
        if all(sum([int(kmer[p] != other[p]) for p in range(k)]) >= 2 for other in kmer_list if kmer != other):
            strong_kmers.append(kmer)
        else:
            weak_kmers.append(kmer)

    return strong_kmers, weak_kmers

# print(pairwise_hamm_dist(['ACAGATTC','ACAGATCC','AGAGATTC']))

def extract_kmers(genome, k):
    unique_kmers = []
    for i in range(len(genome) - k + 1):
        kmer = genome[i:i+k]
        unique_kmers.append(kmer)
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
    output = output+kmer_list
    output = sorted(output)
    # print(len(output))
    return output

def rev_comp(kmer):
    nt_comp = {'A': 'T', 
               'C': 'G',
               'G': 'C',
               'T': 'A'}
    rev_com = ''.join([nt_comp[i] for i in kmer])[::-1]
    return rev_com

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

# ukmers = []
# ukmers += extract_kmers(sample_genomes['Genome_6'], 5)
# su_kmers = quarter_suk(ukmers)
# print(sorted(su_kmers))
    
def run_full_alg_quarter(input, k):
    su_kmers = {}
    genome_kmers = []
    i = 1
    for key in input:
        # print(i,key)
        genome_kmers += extract_kmers(input[key], k)
        print(len(genome_kmers))

    su_kmers = quarter_suk(genome_kmers)
    
    # kmer_freq = Counter(kmer for key in str_kmers for kmer in str_kmers[key])
    # for key in str_kmers:
    #     su_kmers[key] = [kmer for kmer in str_kmers[key] if kmer_freq[kmer] == 1]
    # for key in su_kmers:
    #     print(f'{key} : {len(su_kmers[key])}')
    print('suks' , len(su_kmers))
    return su_kmers

def run_suk(i):
    ukmers = extract_kmers(sample_genomes['Genome_'+str(i+1)], 5)
    full_kmer_list = get_rev_com_kmers(ukmers)
    skmers = quarter_suk(full_kmer_list)
    
    # print(sample_genomes['Genome_'+str(i+1)])
    print(f'Genome_{i+1} :', len(skmers))
    print(sorted(skmers))
    return skmers

# output = []
# for i in range(10):
#     output += [len(run_suk(i))]
# print(output)

# print('Genome: ',run_suk(1))

# start = time.time()
# strongly_unique_kmers = run_full_alg_quarter(sample_genomes,3)
# end = time.time()
# print(end-start)

# times = []
# for i in [50,100,150,200,250]:
#     sample_genomes = {}
#     for record in SeqIO.parse(f'Datasets\diff_sim_datasets\sample_genomes_5k_{i}.fasta', 'fasta'):
#             sample_genomes[record.id] = str(record.seq)
    
#     start = time.time()
#     strongly_unique_kmers = run_full_alg_quarter(sample_genomes,31)
#     end = time.time()
#     print(end-start)
#     times.append(end-start)

# print(times)

plt.plot([15,16,17,18,19,20,21,22,23,24], [26, 30, 34, 37, 40, 43, 46, 49, 52, 55], '.', linestyle = '-')
plt.plot([15,16,17,18,19,20,21,22,23,24], [5, 6, 7, 8, 9, 10, 11, 12, 13, 14], '.', linestyle = '-')
plt.plot([15,16,17,18,19,20,21,22,23,24], [41, 47, 52, 58, 66, 73, 80, 87, 94, 101], '.', linestyle = '-')
plt.plot([15,16,17,18,19,20,21,22,23,24], [0, 0, 0, 1, 2, 3, 4, 5, 6, 7], '.', linestyle = '-')
plt.plot([15,16,17,18,19,20,21,22,23,24], [57, 64, 71, 78, 86, 93, 100, 108, 115, 122], '.', linestyle = '-')
plt.title("Number of SUKs with Varying K in Mycoplasm Genitalium")
plt.xlabel("Value of k")
plt.xticks([15,16,17,18,19,20,21,22,23,24])
plt.ylabel("Number of strongly unique kmers")
plt.legend(['Genome ' + str(p) for p in range(1,6)])
plt.savefig('number_of_suks_new.svg', format = 'svg')
plt.close()
# plt.show()


plt.plot([50,100,150,200,250], [54.835920333862305, 131.00817394256592, 197.56739115715027, 219.81739044189453, 253.0976438522339],marker='o', linestyle='-')
plt.plot([50,100,150,200,250], [27.158747911453247, 67.99119591712952, 108.55513334274292, 168.58467626571655, 199.24515199661255], marker='o', linestyle='-')
plt.plot([50,100,150,200,250], [1.6231374740600586, 3.472759246826172, 5.286814212799072, 7.315382719039917, 7.103714942932129], marker='o', linestyle='-')

plt.title("Runtime v/s Pangenome Similarity (Size: 2kbp)")
plt.xlabel("Number of Mutations")
plt.ylabel("Runtime(s)")
plt.xticks([50,100,150,200,250])
plt.legend(['Alg 1','Alg 2','Alg 3'])
plt.savefig('Runtime_varying_sim.svg', format = 'svg')
# Show the plot
plt.show()

plt.plot([1,2,5,10,15], [5.19626259803772,19.602827787399292,137.64006185531616, 778.8604757785797, 3699.55049514770524], '.', linestyle = '-')
plt.plot([1,2,5,10,15], [2.7661685943603516, 8.387086153030396, 50.419358253479004, 212.74148654937744, 472.36549139022827], '.', linestyle = '-')
plt.plot([1,2,5,10,15], [0.2210984230041504, 0.6606917381286621, 3.4074020385742188, 17.298425912857056, 54.85420346260071], '.', linestyle = '-')
plt.title("Runtime with Varying Genome Size")
plt.xlabel("Genome size (kbp)")
plt.xticks([1,2,5,10,15])
plt.ylabel("Runtime (s)")
plt.legend(['Alg 1', 'Alg 2', 'Alg 3'])
plt.savefig('Runtime_varying_size_new.svg', format = 'svg')
plt.close()





