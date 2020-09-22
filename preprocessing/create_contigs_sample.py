from Bio import SeqIO
import numpy as np


def create_sample_of_100(file, source=False):
    if source:
       fasta_sequences = SeqIO.parse(open(file),'fastq')
    else:
        fasta_sequences = SeqIO.parse(open(file),'fasta')
    counter = 0
    short_samples = []
    for record in fasta_sequences:
        if counter > 100:
            break
        # sample only shorter then 500bp
        if len(record.seq) < 500:
            short_samples.append(record)
            # count  the number of contigs
            counter = counter + 1
    SeqIO.write(short_samples, "../main-code/sample_contigs_k24.contigs.fa", "fasta")

create_sample_of_100("../../contigs-outputs/basic_k-mer24/basic_try_k-mer24.contigs.fa")
