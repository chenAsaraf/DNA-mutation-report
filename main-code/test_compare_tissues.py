from compare_tissues import compare_tissues
import sys


def run_compare_tissues():
    """run compare_tissues
    
    used:
    python3.6 compare_tissues helathy_file_path tumor_file_path output_prefix (optional) test_num (optional)
    helathy_file_path, tumor_file_path - contigs file in FASTA format
    test_num 
    
    """

    healthy = "../../contigs-outputs/healthy/basic_k-mer24/basic_try_k-mer24.contigs.fa"
    tumor = "../../contigs-outputs/tumor/basic_k-mer24_T/basic_k-mer24_T.contigs.fa"
    argvlen = len(sys.argv)
    if argvlen < 3:  # missing contigs file path
        # raise error
    else:
        healthy = sys.argv[1]
        tumor = sys.argv[2]
    if argvlen > 3:
        output_prefix = sys.argv[3]
    else: 
        output_prefix = "output"
    if argvlen > 4:
        test_num = sys.argv[4]
    else:
        test_num = None

    compare_tissues(healthy, tumor, output_prefix, test=False, test_num=test_num)

test_compare_tissues()
