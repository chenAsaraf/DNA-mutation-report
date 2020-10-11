# DNAMutationReport
##### A Tool for analysis of tumor mutation burden

We present here the program and the results of our research for constructing a feasibly running algorithm that once given a collection of unorganised reads (strings of nucleobases A, C, G, T) from a tumour tissue and a healthy tissue of the same specimen, the algorithm searches for mutations of certain types. Which, in turn, when evaluating the number of mutations, categorizing mutations and identifying precise changes in comparison to a normal genome, can lead to better detection, diagnosis and treatment provided healthcare.
The algorithm uses the reads to partially assemble the genome, with the help of an external existing de-novo assembly program. Once the reads are mapped into longer sequences, called ‘contigs’, the algorithm uses a dictionary type data-structure to reduce the number of comparisons between them.

## Tools we used:
* **Minia assembler** - A short-read assembler based on a de Bruijn graph, the output is a set of contigs. see more at [Minia page](https://github.com/GATB/minia).

minia command line we used:
```
./minia -in reads.fa -kmer-size 24 -abundance-min 3 -out output_prefix
```
The main parameters are:
  1. reads.fa* – the input file(s)
  2. kmer-size 24 – k-mer length (integer), the number may vary depending on user choice
  3. abundance-min 3 - hard cut-off to remove likely erroneous, low-abundance k-mers
  4. output_prefix – any prefix string to store output contigs as well as temporary files for this assembly

* **edit-distance** - Python module for computing edit distances and alignments between sequences. see more at [edit-distancw page](https://github.com/belambert/edit-distance).

## How to Use:

in 'main-code' folder run:
```
python3.6 run_compare_tissues helathy_file_path tumor_file_path output_prefix(optional) test(optional) test_num(optional)
```
The parameters are:
*helathy_file_path, tumor_file_path - contigs file in FASTA format
*test - this variable designed to assist in the software
         development process. If 'test' argument exist then the
         software will only run up to test_num contigs.
*test_num - int (optional) this parameter used only in case 'test' argument 
               exist
