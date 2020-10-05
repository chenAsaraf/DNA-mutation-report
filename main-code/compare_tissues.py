from collections import defaultdict
from Bio import SeqIO
from build_dictionary import tissueDictionary
from matplotlib import pyplot as plt
from mutations_distance import PointMutation
from contigs_analysis import filter_contigs_by_size
import time
# 1) parse the tumor_cell with non-overlaping window and find there bucket in the dictionary
# 2) search for each contig in the bucket-list in the BST
# 3) find alignment (maby another function)
# 4) send to another function that find the mutation by edit distance

"""
"""
def find_similar_section(tumor_file, output_prefix, k, dictionary, healthyStorage, test=False, test_num=1000):
    # for searching the correct bucket in the dictionary
    # run throw all NON-overlaping windows of length k in the sequence (all k-mer)
    filtered_tumor_file, num_tumor_contigs = filter_contigs_by_size(tumor_file, 'filtered_tumor_contigs', test=test, test_num=test_num)
    print("number of filtered tumor contigs:", num_tumor_contigs)
    # initialize the object to save the mutations
    mutations_report = PointMutation(output_prefix)
    records = SeqIO.parse(open(filtered_tumor_file), 'fasta')
    for tumor_seq in records:
        contig_len = len(str(tumor_seq.seq))
        for window in range(0, contig_len - k, k): # (last k argument - to jump k chars each step)
            if str(tumor_seq.seq)[window: window + k] in dictionary.keys():
                # Go through all the strings of the healthy tissue in this bucket
                healthy_tissue_records = dictionary[str(tumor_seq.seq)[window: window + k]]
                for record in healthy_tissue_records:
                    healthy_seq = healthyStorage[record.id]
                    # For each alignment - find the overlapping parts and send to the Edit-Distance function
                    for healthy_idx in record.indexes:
                        healthy, tumor = find_overlap(str(healthy_seq), str(tumor_seq.seq), healthy_idx, window)
                        mutations_report.editDistance(tumor, healthy)
    return mutations_report

def find_overlap(healthy_seq, tumor_seq, healthy_idx, tumor_idx):
    begin_healthy = end_healthy = begin_tumor = end_tumor = 0
    # 1) Find the beginning of the overlap:
    # If one or both of the indexes is 0 - the start of the sequence, leave the indexes the same
    if tumor_idx == 0 or healthy_idx == 0:
        begin_healthy = healthy_idx
        begin_tumor = tumor_idx
    # Else- if one is larger then the other - take the sequence with the smaller index
    # from the beginning, and take the sequence with the larger index from:
    # larger_index - smaller_index
    else:
        if healthy_idx < tumor_idx:
            begin_healthy = 0
            begin_tumor = tumor_idx - healthy_idx
        else:
            begin_tumor = 0
            begin_healthy = healthy_idx - tumor_idx
    # 2) Find the ending of the overlap:
    remaining_healthy = len(healthy_seq) - healthy_idx
    remaining_tumor = len(tumor_seq) - tumor_idx
    if remaining_healthy < remaining_tumor:
        end_healthy = healthy_idx + remaining_healthy
        end_tumor = tumor_idx + remaining_healthy
    else:
        end_healthy = healthy_idx + remaining_tumor
        end_tumor = tumor_idx + remaining_tumor
    return healthy_seq[begin_healthy : end_healthy], tumor_seq[begin_tumor : end_tumor]


def compare_tissues(healthy_file, tumor_file, output_prefix, test=False, test_num=1000):
    start_time_build_dict = time.time()

    if test:
        dictBuilder = tissueDictionary(healthy_file, test=True, test_num=test_num)
    else:
        dictBuilder = tissueDictionary(healthy_file)

    print("--- %s seconds to BUILD dictionary ---" % (time.time() - start_time_build_dict))

    dictionary, contigsStorage = dictBuilder.get_dictionary_and_storage()
    k = dictBuilder.getK()

    start_time_compare = time.time()

    mutations_report = find_similar_section(tumor_file, output_prefix, k, dictionary, contigsStorage, test=test, test_num=test_num)

    print("--- %s seconds to COMPARE all tissues ---" % (time.time() - start_time_compare))

    print()
    print("MUTATIONS REPORT:")
    print("INSERTS:", mutations_report.inserts)
    print("REPLACES:", mutations_report.replaces)
    print("DELETES:", mutations_report.deletes)
    print("NUMBER OF CONTIGS WHERE COMPARED:", mutations_report.counterOfCompares)
    avg_length_of_contig = mutations_report.sumOfLength/mutations_report.counterOfCompares
    print("AVG LENGTH OF CONTIG:", avg_length_of_contig)
    print("AVG NUMBER OF INSERTS:", mutations_report.counters[0]/mutations_report.counterOfCompares)
    print("AVG NUMBER OF REPLACES:", mutations_report.counters[1]/mutations_report.counterOfCompares)
    print("AVG NUMBER OF DELETES:", mutations_report.counters[2]/mutations_report.counterOfCompares)
    # Creating plot
    fig = plt.figure(figsize=(10, 7))
    plt.pie(mutations_report.counters, labels=["inserts", "replaces", "deletes", "matches"], autopct='%1.1f%%')

    # save plot
    fig.savefig(output_prefix + ".png")
