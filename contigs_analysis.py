from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt


def plot_histogram(array_of_lengthes):
    bins = np.concatenate(([100,200,300,400,500,1000],np.arange(2000,6000,1000)), axis=None)
    fig = plt.figure()
    plt.hist(array_of_lengthes, bins=bins, density=True, facecolor='g')
    plt.title("Histogram of contgis lengthes - with k-mer : 24 for R2")
    plt.xlabel('Length (bp)')
    plt.grid(True)
    fig.savefig('histogram_basic24_R2.png')


def read_fasta(contigs_file):
    fasta_sequences = SeqIO.parse(open(contigs_file),'fasta')
    counter = 0
    min_length = float('inf')
    max_length = float('-inf')
    sum_length = 0
    list_of_length =[]
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        # count  the number of contigs
        counter = counter + 1
        contig_len = len(sequence)
        list_of_length.append(contig_len)
        # sum the length of the contigs:
        sum_length = sum_length + contig_len
        # find the minimum and the maximum length of contigs:
        if contig_len < min_length:
            min_length = contig_len
        if contig_len > max_length:
            max_length = contig_len
    # find the avarege length of contigs
    avg = sum_length / counter
    """print("number of contigs:", counter)
    print("min length of contig:", min_length)
    print("max length of contig:", max_length)
    print("avarege length of contigs:", avg)"""
    plot_histogram(np.array(list_of_length))
    """print("histogram:")
    print(hist)
    print("bin_edges:")
    print(bin_edges)"""

def read_source_file(reads_file):
	count = 0
	for record in SeqIO.parse(open(reads_file), "fastq-illumina"):
		count = count + 1
	print("the number of reads are : ", count)

read_fasta("../contigs-outputs/basic_k-mer24-R2/basic_k-mer24-R2.contigs.fa")
read_source_file("../source_files+minia/sample_TB0001955-16933-N_R1_001.fastq") 
