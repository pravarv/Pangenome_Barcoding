# Pangenome_Barcoding

The python files contain implementations of the code for the 3 algorithms mentioned in the project report, along with an implementation of the Quarter Alg from Rahmann et al. (2024).

Implementation of Algorithm 1 (Brute Force) : alg_1_genome_wise.py

Implementation of Algorithm 2 (with PanKmer): alg_2_wPankmer.py

Implementation of Algorithm 3 (Pankmer + Quarter): alg_3_pankmer_quarter.py

Implementation of Quarter Algorithm : quarter_code.py

Self Dataset Generator : pangenome_generator.py

The zip files are the datasets used for debugging and generating relevant results. Biological Data.zip contains 3 pangenomes, obtained from https://www.ncbi.nlm.nih.gov/datasets/genome/. 

Mycoplasm Genitalium (5 genomes, 580 kbp)

Brucella Melitensis (5 genomes, 3.3 Mbp)

Sacchromyces Cerevisiae (5 genomes, 12Mbp)

The self-generated datasets include debugging and testing datasets. 

'_diff_sim_datasets_' contains pangenomes with same genome size but varying similarity.

'_diff_size_datasets_' contains pangenomes with different genome size but a constant degree of similarity.

The remaining fasta files are sample pangenomes with file name indicating the properties of the pangenome. Filename : 'sample_genomes_{_number_}_{_genome size_}_s{_number of mutations_}.fasta



