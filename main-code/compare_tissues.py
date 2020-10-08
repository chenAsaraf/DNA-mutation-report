from collections import defaultdict
from Bio import SeqIO
from matplotlib import pyplot as plt
import time
import numpy as np
from build_dictionary import TissueDictionaryBuilder
from mutations_distance import PointMutation
from contigs_analysis import filter_contigs_by_size
import plot_diagrams as diagrams


# 1) parse the tumor_cell with non-overlaping window and find there bucket in the dictionary
# 2) search for each contig in the bucket-list in the BST
# 3) find alignment (maby another function)
# 4) send to another function that find the mutation by edit distance

"""
"""
def find_similar_section(tumor_file, output_prefix, k, dictionary, healthyStorage, test=False, test_num=1000):
    # for searching the correct bucket in the dictionary
    # run throw all non-overlaping windows of length k in the sequence (all k-mer)
    filtered_tumor_file, num_tumor_contigs = filter_contigs_by_size(tumor_file, 'filtered_tumor_contigs', test=test, test_num=test_num)
    print("number of filtered tumor contigs:", num_tumor_contigs)
    # initialize the object to save the mutations
    mutations_report = PointMutation(output_prefix)
    records = SeqIO.parse(open(filtered_tumor_file), 'fasta')
    statistics = []
    for tumor_seq in records:
        statistics_of_compares = 0
        contig_len = len(str(tumor_seq.seq))
        for window in range(0, contig_len - k, k): # (last k argument - to jump k chars each step)
            if str(tumor_seq.seq)[window: window + k] in dictionary.keys():
                # Go through all the strings of the healthy tissue in this bucket
                healthy_tissue_records = dictionary[str(tumor_seq.seq)[window: window + k]]
                for record in healthy_tissue_records:
                    healthy_seq = healthyStorage[record.id]
                    # For each alignment - find the overlapping parts and send to the Edit-Distance function
                    for healthy_idx in record.indexes:
                        # Increase the number of compared contigs (For later analysing)
                        statistics_of_compares = statistics_of_compares + 1
                        healthy, tumor = find_overlap(str(healthy_seq), str(tumor_seq.seq), healthy_idx, window)
                        mutations_report.editDistance(tumor, healthy)
        statistics.append(statistics_of_compares)

    return mutations_report, statistics

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


def calc_mutations_amount(mutations_report):
    """ Sum all the characters that changed per mutation type """

    inserts_amount = 0
    for char, number in mutations_report.inserts.items():
        inserts_amount = inserts_amount + number
    replaces_amount = 0
    for chars, number in mutations_report.replaces.items():
        replaces_amount = replaces_amount + number
    deletes_amount = 0
    for char, number in mutations_report.deletes.items():
        deletes_amount = deletes_amount + number
    return inserts_amount, replaces_amount, deletes_amount


def calculate_percentages(mutations_report, inserts_amount, replaces_amount, deletes_amount):
    """ Calculate for each mutation type the percentage of characters that changed
    out of all the characters of the diseased tissue we tested"""

    mutations_report.set_in_percentages(
        (float(inserts_amount) / float(mutations_report.sumOfLengths)) * 100 )
    mutations_report.set_rep_percentages(
        (float(replaces_amount) / float(mutations_report.sumOfLengths)) * 100 )
    mutations_report.set_del_percentages(
        (float(deletes_amount) / float(mutations_report.sumOfLengths)) * 100 )


def compare_tissues(healthy_file, tumor_file, output_prefix, test=False, test_num=1000):
    """ Main Function """

    # To measure program times:
    start_time_build_dict = time.time()
    # Build the dictionary from the healthy tissue file
    if test:  # While writing the program we will run with the 'test'=True parameter
        dict_builder = TissueDictionaryBuilder(healthy_file, test=True, test_num=test_num)
    else:
        dict_builder = TissueDictionaryBuilder(healthy_file)
    # Prints to stdout for indication of program progress
    print()
    print("---------------------------------------")
    print("%.2f seconds to BUILD dictionary " % (time.time() - start_time_build_dict))
    print("---------------------------------------")
    print()
    # Get the dictionary and the contigs storage
    dictionary, contigs_storage = dict_builder.get_dictionary_and_storage()
    k = dict_builder.get_k()  # The dictionary window size
    start_time_compare = time.time()
    print()
    print("start to compare tissues...")
    # Compare the healthy tissue to the tumor tissue, and return the report
    mutations_report, statistics = find_similar_section(tumor_file, output_prefix, k, dictionary, contigs_storage, test=test, test_num=test_num)
    print()
    print("---------------------------------------")
    print("%.2f seconds to COMPARE all tissues " % (time.time() - start_time_compare))
    print("---------------------------------------")
    print()
    # Calculate statistics of comparisons:
    sum_comparisons = sum(statistics)
    avg_compares_statistics = sum_comparisons/len(statistics)
    print("number of all comparison:", sum_comparisons)
    print("statistic of contigs comparison per one: ", avg_compares_statistics)
    # Calculate statistics per mutation type:
    inserts_amount, replaces_amount, deletes_amount = calc_mutations_amount(mutations_report)
    calculate_percentages(mutations_report, inserts_amount, replaces_amount, deletes_amount)
    # Print results to the screen:
    print(" ~ Mutations Report ~ ")
    print()
    print("The cancerous tissue from: ", tumor_file)
    print("The noraml tissue from:", healthy_file)
    print()
    print("~ Inserts Amount:", inserts_amount,
          ", Percentage of all characters: %.2f " % mutations_report.in_percentages, "%")
    print("\t \t", mutations_report.inserts)
    print("~ Replaces Amount:", replaces_amount,
          ", Percentage of all characters: %.2f " % mutations_report.rep_percentages, "%")
    print("\t \t", mutations_report.replaces)
    print("~ Deletes Amount:", deletes_amount,
          ", Percentage of all characters: %.2f" % mutations_report.del_percentages, "%")
    print("\t \t", mutations_report.deletes)
    print("~ Number of comparisons actually entered into the report:",
          mutations_report.counterOfComparisons)
    # Create diagrams:
    diagrams.general_diagram(inserts_amount, replaces_amount, deletes_amount, output_prefix)
    diagrams.inserts_diagram(mutations_report, output_prefix)
    diagrams.replaces_diagram(mutations_report, output_prefix)
    diagrams.deletes_diagram(mutations_report, output_prefix)
    # Create Histogram of all distances for optimizations
    title = "Histogram of all comparisons distances"
    diagrams.distance_histogram(np.array(mutations_report.listOfDistances), "distances-histogram-" + output_prefix, title)

