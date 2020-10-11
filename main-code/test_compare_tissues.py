from compare_tissues import compare_tissues
import sys


def run_compare_tissues():
    """run compare_tissues
    
    used:
    python3.6 run_compare_tissues helathy_file_path tumor_file_path output_prefix (optional) test (optional) test_num (optional)
    helathy_file_path, tumor_file_path - contigs file in FASTA format
    test - this variable designed to assist in the software
             development process. If 'test' argument exist then the
             software will only run up to test_num contigs.
    test_num - int (optional) this parameter used only in case 'test' argument 
               exist
    
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
        if sys.argv[4] == 'test':
            test = True
            test_num = sys.argv[5]
    else:
        test = False
        test_num = None

    compare_tissues(healthy, tumor, output_prefix, test=test, test_num=test_num)

run_compare_tissues()
