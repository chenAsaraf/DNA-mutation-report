from Bio import SeqIO
import time
import numpy as np
from collections import defaultdict
from matplotlib import pyplot as plt
from build_dictionary import TissueDictionaryBuilder
from mutations_distance import PointMutation
from contigs_analysis import filter_contigs_by_size
import plot_diagrams as diagrams

"""Main file contains the main functions to run the program

The main function is 'compare_tissues'.
Program results are:
    - Diagrams: for each mutations type and a distances histogram for 
            analysing in developing process
    - 'txt' file of sample strings from the contigs with farther details
    - PointMutation object, contains the results, and print to stdout 
            its fields
            
            
The program stages:
    1) Parse the healthy tissue contigs and store their pointer in the
        dictionary with TissueDictionaryBuilder object
    2) For each contig from the tumor tissue contigs:
        parse the sequence with non-overlapping window and find their entries
        in the dictionary (find_similar_section function)
        3) For each contig pointer in the each entry: 
            4) Find the alignment between the tumor sequence to the healthy
                sequence, and cut only the overlapping part (find_overlap)
            5) Calculate the Edit-Distance between the overlapping part 
                with PointMutation object
"""


def find_similar_section(tumor_file, output_prefix, k, dictionary, healthy_storage,
                         test=False, test_num=1000):
    """ For each contig in the tumor file move the window along the sequence and
        find the similar contig from the healthy file. then compare their
        overlapping parts with edit distance

        :param tumor_file: 'FASTA' file
                    tumor tissue contigs file
        :param output_prefix: string
                    prefix of the output files
        :param k: int
                    the window size (k-mer) of the dictionary-builder
        :param dictionary:  defaultdict(list)
                    holds the k-mer as keys and list of dictionaryItem
                    as values
        :param healthy_storage: list
                     holds the healthy tissue contigs sequence in
                     the index of their id
        :param test: boolean (optional, default = False)
                    a variable designed to assist in the software
                    development process. If test=True then the
                    software will only run up to test_num contigs.
        :param test_num: int (optional, default = 1000)
                    this parameter used only in case test=True
        :return mutations_report: PointMutation object
        :return statistics: list
                    list of comparisons number for each tumor contig
    """
    # Select the contigs from the tumor tissue in a certain length range
    filtered_tumor_file, num_tumor_contigs = filter_contigs_by_size(
        tumor_file, 'filtered_tumor_contigs', test=test, test_num=test_num)
    print("number of filtered tumor contigs:", num_tumor_contigs)
    # Initialize object to save the mutations
    mutations_report = PointMutation(output_prefix)
    # Open the filtered tumor contigs file
    records = SeqIO.parse(open(filtered_tumor_file), 'fasta')
    # Statistic of comparisons for each sequence in tumor tissue:
    statistics = []
    for tumor_seq in records:
        statistics_of_compares = 0
        contig_len = len(str(tumor_seq.seq))
        # Move the window along the contig sequence
        for window in range(0, contig_len - k, k):
            if str(tumor_seq.seq)[window: window + k] in dictionary.keys():
                healthy_tissue_records = dictionary[str(tumor_seq.seq)[window: window + k]]
                # Go through the list of all contigs of the healthy tissue in this bucket
                for record in healthy_tissue_records:
                    healthy_seq = healthy_storage[record.id]
                    # For each index of the window in the healthy contig -
                    for healthy_idx in record.indexes:
                        # Increase the number of compared contigs (For later analysing)
                        statistics_of_compares = statistics_of_compares + 1
                        # Find the overlapping parts
                        healthy, tumor = find_overlap(str(healthy_seq), str(tumor_seq.seq), healthy_idx, window)
                        # Send to the Edit-Distance function
                        mutations_report.editDistance(tumor, healthy)
        statistics.append(statistics_of_compares)
    return mutations_report, statistics


def find_overlap(healthy_seq, tumor_seq, healthy_idx, tumor_idx):
    """Find the overlapping parts of the contigs sequences

    :param healthy_seq: string
                the contig sequence from the healthy tissue
    :param tumor_seq: string
                the contig sequence from the tumor tissue
    :param healthy_idx: int
                the starting index of the window in the healthy tissue
    :param tumor_idx:
                the starting index of the window in the tumor tissue
    :return: 2 strings in the same length from the both tissues
                represents the overlapping parts
    """
    # 1) Find the beginning of the overlap:
    #   if one or both of the indexes is 0 - the start of the sequence,
    #   leave the indexes the same
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
    return healthy_seq[begin_healthy: end_healthy], tumor_seq[begin_tumor: end_tumor]


def calc_mutations_amount(mutations_report):
    """Sum all the characters that changed per mutation type

    :param mutations_report: PointMutation object
    :return inserts_amount: int
                total inserts in the report
    :return replaces_amount: int
                total replaces in the report
    :return deletes_amount:
                total deletes in the report
    """
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
    """Calculate for each mutation type the percentage of characters that changed
    out of all the characters of the diseased tissue we tested

    :param mutations_report: PointMutation object
    :param inserts_amount: int
                total inserts in the report
    :param replaces_amount:int
                total replaces in the report
    :param deletes_amount:int
                total deletes in the report
    """
    mutations_report.set_in_percentages(
        (float(inserts_amount) / float(mutations_report.sumOfLengths)) * 100)
    mutations_report.set_rep_percentages(
        (float(replaces_amount) / float(mutations_report.sumOfLengths)) * 100)
    mutations_report.set_del_percentages(
        (float(deletes_amount) / float(mutations_report.sumOfLengths)) * 100)


def compare_tissues(healthy_file, tumor_file, output_prefix, test=False, test_num=1000):
    """  Main Function - generate diagrams for each mutations type,
    a 'txt' file of sample strings from the contigs with farther details
    and print to stdout the results

    :param healthy_file: 'FASTA' file
                healthy tissue contigs file
    :param tumor_file: 'FASTA' file
                tumor tissue contigs file
    :param output_prefix: string
                prefix of the output files
    :param test: boolean (optional, default = False)
                a variable designed to assist in the software
                development process. If test=True then the
                software will only run up to test_num contigs.
    :param test_num: int (optional, default = 1000)
                this parameter used only in case test=True
    :return: mutations_report: PointMutation object
                contains the final results
    """

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
    return mutations_report
