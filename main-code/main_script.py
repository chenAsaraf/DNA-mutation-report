from collections import defaultdict
from Bio import SeqIO
import build_dictionary

helthy_cell = ""
tumor_cell=""
dictionary = build_dictionary(helthy_cell)
serach_alignment(tumor_cell, dictionary)


