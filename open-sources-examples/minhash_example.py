from snapy import MinHash, LSH
import numpy as np
from fasta_parser import parse_to_list

contigs_file = "../contigs-outputs/basic_k-mer24/basic_try_k-mer24.contigs.fa"
print("parsing contigs to list...")
contigs_list, num_contigs = parse_to_list(contigs_file)
print("number of contigs (shorter then 1,00bp):", num_contigs)
labels = np.arange(num_contigs).tolist()
# Create MinHash object.
print("creating minhash object...")
minhash = MinHash(contigs_list, n_gram=24)
# Create LSH model.
print("creating LSH object...")
lsh = LSH(minhash, labels)
print("query object:")
print("similar strings to the string in index 0:", lsh.query(0))
print("similar strings to the string in index 1:", lsh.query(1))
print("similar strings to the string in index 2:", lsh.query(2))
print("similar strings to the string in index 3:", lsh.query(3))
print("similar strings to the string in index 4:", lsh.query(4))
print("similar strings to the string in index 5:", lsh.query(5))
print(lsh.contains())

