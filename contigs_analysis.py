from Bio import SeqIO


def some_function(name, sequence):
    print(name, sequence)

contigs_file = "../contigs-outputs/basic_k-mer24/basic_try_k-mer24.contigs.fa"
fasta_sequences = SeqIO.parse(open(contigs_file),'fasta')
counter = 0
min_length = float('inf')
max_length = float('-inf')
sum_length = 0
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    #new_sequence = some_function(name, sequence)
    # count  the number of contigs
    counter = counter + 1
    contig_len = len(sequence)
    # sum the length of the contigs:
    sum_length = sum_length + contig_len
    # find the minimum and the maximum length of contigs:
    if contig_len < min_length:
        min_length = contig_len
    if contig_len > max_length:
        max_length = contig_len

# find the avarege length of contigs
avg = sum_length / counter

print("number of contigs:", counter)
print("min length of contig:", min_length)
print("max length of contig:", max_length)
print("avarege length of contigs:", avg)
