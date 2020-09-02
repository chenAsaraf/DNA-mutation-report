from snapy import MinHash, LSH
import numpy as np
from fasta_parser import parse_to_list

contigs_file = "../contigs-outputs/basic_k-mer24/basic_try_k-mer24.contigs.fa"
contigs_list = parse_to_list(contigs_file)
#contigs_list = ['AAAACGCAAAACGC', 'AAAACTCAAAACTC', 'AAACTCTAAACTCT',
#                'CCTCTCCCTCTC', 'GTCAAGTGGTCAAGTG', 'GTCATTTGTCATTT']
labels = np.arange(len(contigs_list)).tolist()
print(labels)
# Create MinHash object.
minhash = MinHash(contigs_list, n_gram=9)
# Create LSH model.
lsh = LSH(minhash, labels)
print("similar strings to the string in index 0:", lsh.query(0))
print("similar strings to the string in index 1:", lsh.query(1))
print("similar strings to the string in index 2:", lsh.query(2))
print("similar strings to the string in index 3:", lsh.query(3))
print("similar strings to the string in index 4:", lsh.query(4))
print("similar strings to the string in index 5:", lsh.query(5))
print(lsh.contains())

