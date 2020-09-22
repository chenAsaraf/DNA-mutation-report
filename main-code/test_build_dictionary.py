from collections import defaultdict
from Bio import SeqIO
from build_dictionary import cellDictionary


helthy_cell = "sample_contigs_k24.contigs.fa"
dictionary, contigsStorage = cellDictionary(helthy_cell).get_dictionary_and_storage()
print()
print("test:")
for kmer, contigs in dictionary.items():
    print(kmer, "->")
    for contigData in contigs:
        print("\t", contigData.id)
        for index in contigData.indexes:
            print("\t \t", index)

