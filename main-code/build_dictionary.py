from collections import defaultdict
from Bio import SeqIO

"""
Build the dictionary from the healthy cell, for later use
"""

class build_dictionary:

    def __init__(self, contigs_file):
        self.dictionary = defaultdict(list)
        records = SeqIO.parse(open(contigs_file), 'fasta')
        for contig in records:
            # TODO: save the contig in BST
            # parse the windows to dictionary:
            self.parse_window(contig)


    def parse_window(self, contig):
    # run throw all overlaping windows of length 10 in the sequence (all 10-mer)
        for i in range(len(str(contig.seq))-10):
            # key: kmer
            # value: pairs list of (contig id, index of window)
            self.dictionary[str(contig.seq)[i : i+10]].append((contig.id, i))
