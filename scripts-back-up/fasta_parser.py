from Bio import SeqIO


def parse_to_list(contigs_file):
    fasta_sequences = SeqIO.parse(open(contigs_file),'fasta')
    contigs_list = []
    counter = 0
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        if len(sequence) < 1000:
            contigs_list.append(sequence)
            counter = counter + 1
    return contigs_list, counter
