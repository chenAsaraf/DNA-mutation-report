from collections import defaultdict
from Bio import SeqIO
import time
from build_dictionary import TissueDictionaryBuilder

helthy_cell = "../../contigs-outputs/basic_k-mer24/basic_try_k-mer24.contigs.fa"
start_time_build_dict = time.time()
dictionary, contigsStorage = TissueDictionaryBuilder(helthy_cell, test=True).get_dictionary_and_storage()
print("--- %s seconds to build dictionary ---" % (time.time() - start_time_build_dict))
print()
print("Analysis of dictionary bucket's size:")
lengths_list = []
min_length = float('inf')
max_length = float('-inf')
start_time_iterate_dict = time.time()
for kmer, contigs in dictionary.items():
    # print(kmer, "->")
    counter = 0
    for contigData in contigs:
        counter = counter + 1
        # print("\tcontigs:")
        # print("\t-", contigData.id)
        # print("\tindexes of contig", contigData.id)
        # for index in contigData.indexes:
            # print("\t\t", index)
    if min_length > counter:
        min_length = counter
    if max_length < counter:
        max_length = counter
    lengths_list.append(counter)

print("--- %s seconds to iterate through all the dictionary ---" % (time.time() - start_time_iterate_dict))
print("minimum number in bucket:", min_length)
print("maximum number in bucket:", max_length)
print("average number in bucket:", (sum(lengths_list) / len(lengths_list)))
