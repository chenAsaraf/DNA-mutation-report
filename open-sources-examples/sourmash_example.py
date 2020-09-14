import sourmash
from Bio import SeqIO

#contigs_file = "../contigs-outputs/basic_k-mer24/basic_try_k-mer24.contigs.fa"
#contigs_file = "contig_sample_k24.contigs.fa"
reads_file = "../source_files+minia/sample_TB0001955-16933-N_R1_001.fastq"
fasta_sequences = SeqIO.parse(open(reads_file),'fasta')
minhashes = []
mh = sourmash.MinHash(n=0, ksize=24, scaled=1)
for contig in fasta_sequences:
    name, sequence = contig.id, str(contig.seq)
    mh.add_sequence(sequence, True)
    minhashes.append(mh)

print()
print()
for i, mh1 in enumerate(minhashes):
    for j, mh2 in enumerate(minhashes):
        sim = mh1.similarity(mh2)
        print("similarity of", i, j, ":" ,sim)

