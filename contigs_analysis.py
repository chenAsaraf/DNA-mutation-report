from Bio import SeqIO


def some_function(name, sequence):
    print(name, sequence)

contigs_file = "../contigs-outputs/basic_k-mer24/basic_try_k-mer24.contigs.fa"
fasta_sequences = SeqIO.parse(open(contigs_file),'fasta')
counterOfContigs = 0
counterOfBigContigs = 0
counterOfHugeContigs = 0
min_length = float('inf')
max_length = float('-inf')
sum_length = 0
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    #new_sequence = some_function(name, sequence)
    # count  the number of contigs
    counterOfContigs = counterOfContigs + 1
    contig_len = len(sequence)
    # sum the length of the contigs:
    sum_length = sum_length + contig_len
    # find the minimum and the maximum length of contigs:
    if contig_len < min_length:
        min_length = contig_len
    if contig_len > max_length:
        max_length = contig_len
    # count the number of contigs with length of 150 - 300
    if contig_len > 150 and contig_len <= 300:
        counterOfBigContigs += 1
    # count the number of contigs with length of over 300
    if contig_len > 300:
        counterOfHugeContigs += 1

# find the avarege length of contigs
avg = sum_length / counter

print("number of contigs:", counterOfContigs)
print("min length of contig:", min_length)
print("max length of contig:", max_length)
print("avarege length of contigs:", avg)
print("number of contigs with length of 150 - 300:", counterOfBigContigs)
print("number of contigs with length of over 300:", counterOfHugeContigs)
