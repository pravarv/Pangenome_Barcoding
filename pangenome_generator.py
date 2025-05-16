import random

length_of_genome = 2000
num_of_mutations = 50
num_of_genomes = 10
# random.seed(42)

ref_genome = ''.join(random.choices(['A','G','C','T'], k = length_of_genome))
print(len(ref_genome))

def insert_mut(genome, sub_num):
    pos_list = random.sample(range(1,len(genome)+1), sub_num)
    for pos in pos_list:
        nts = [x for x in 'ACGT']
        p = genome[pos-1]
        nts.remove(p)
        r_1 = random.randint(0,2)
        genome = genome[:pos-1] + nts[r_1] + genome[pos:]
    
    return genome

def generate_genomes(reference_genome, num_genomes,mut_num):
    genomes = {}
    for i in range(num_genomes):
        genomes[f'Genome_{i+1}'] = insert_mut(reference_genome,mut_num)  
    return genomes

#Creating a fasta file with pangenome for PanKmer Input
for mut_num in [50,100,150,200,250]:
    genomes = generate_genomes(ref_genome,num_of_genomes, mut_num)
    with open(f'sample_genomes_2k_{mut_num}.fasta', "w") as fasta_file:
        for genome_id, sequence in genomes.items():
            fasta_file.write(f">{genome_id}\n")
            fasta_file.write(f"{sequence}\n")


