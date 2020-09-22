from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import math

def plot_histogram(array_of_lengthes):
    bins = np.concatenate(([50,100,200,300,400,500,1000],np.arange(2000,9000,1000)), axis=None)
    fig = plt.figure()
    plt.hist(array_of_lengthes, bins=bins, density=True, facecolor='g')
    plt.title("Histogram of Normal Cell Contigs (kmer:24)")
    plt.xlabel('Length (bp)')
    plt.grid(True)
    fig.savefig('histogram-normal-contigs-k24_1.png')


def analyse_file(file, source=False):
    if source:
       fasta_sequences = SeqIO.parse(open(file),'fastq')
       print("starting with opening source reads file")
    else:
        fasta_sequences = SeqIO.parse(open(file),'fasta')
        print("starting with opening contigs file")
    counter = 0
    min_length = float('inf')
    max_length = float('-inf')
    sum_length = 0
    list_of_length =[]
    for record in fasta_sequences:
        name, seq = record.id, str(record.seq)
        # count  the number of contigs
        counter = counter + 1
        seq_len = len(seq)
        list_of_length.append(seq_len)
        # sum the length of the contigs:
        sum_length = sum_length + seq_len
        # find the minimum and the maximum length of contigs:
        if seq_len < min_length:
            min_length = seq_len
        if seq_len > max_length:
            max_length = seq_len
    # find the avarege length of contigs
    avg = sum_length / counter
    print("number of sequences:", counter)
    print("min length of sequence:", min_length)
    print("max length of sequence:", max_length)
    print("avarege length of sequence:", avg)
    bins = np.concatenate(([100,200,300,400,500,1000],np.arange(2000,9000,1000)), axis=None)
    hist, bin_edges = np.histogram(np.array(list_of_length), bins=bins)
    print("histogram:")
    print(hist)
    print("bin_edges:")
    print(bin_edges)
    #plot_histogram(np.array(list_of_length))


def find_window_size(contigs_file):
    records = SeqIO.parse(open(contigs_file),'fasta')
    for record in records:
        sequence_len = len(str(record.seq))
        min_length = float('inf')
        # find the minimum length
        if sequence_len < min_length:
            min_length = sequence_len
    num_of_illumina_mistakes = math.ceil(min_length / 10)
    print("num_of_illumina_mistakes", num_of_illumina_mistakes)
    num_of_windows = num_of_illumina_mistakes + 1
    bp_per_window = math.floor(min_length / num_of_windows)
    return bp_per_window

#analyse_file("../contigs-outputs/basic_k-mer24/basic_try_k-mer24.contigs.fa")
#analyse_file("../source_files+minia/sample_TB0001955-16933-N_R1_001.fastq", source=True)
#analyse_file("../main-code/sample_contigs_k24.contigs.fa")
#read_source_file("../source_files+minia/sample_TB0001955-16933-N_R1_001.fastq") 
