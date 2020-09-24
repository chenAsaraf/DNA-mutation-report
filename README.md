# DNAMutationReport
## Tool for analysis of tumor mutation burden

##### Tools we used:
- Minia assembler: 
minia command line we used:
'''
./minia -in reads.fa -kmer-size 24 -abundance-min 3 -out output_prefix
'''
The main parameters are:
1. **reads.fa** – the input file(s)
2. **kmer-size 24 ** – k-mer length (integer), the number may vary depending on user choice
3. **abundance-min 3** - hard cut-off to remove likely erroneous, low-abundance k-mers
4. **output_prefix** – any prefix string to store output contigs as well as temporary files for this assembly
see more at [Minia page](https://github.com/GATB/minia)
