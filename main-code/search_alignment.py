from collections import defaultdict
from Bio import SeqIO
import build_dictionary

# 1) parse the tumor_cell with non-overlaping window and find there bucket in the dictionary
# 2) search for each contig in the bucket-list in the BST
# 3) find alignment (maby another function)
# 4) send to another function that find the mutation by edit distance

