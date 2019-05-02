import pysam
import argparse

parser = argparse.ArgumentParser(description='Utility for converting bam file to fastq file.')
parser.add_argument('input_bam_file', type=str, nargs='?', help='provide input bam file path here')
parser.add_argument('out_file', type=str, nargs='?', help='provide output fastq file path here')
args = parser.parse_args()

outf = open(args.out_file, "w")
samfile = pysam.AlignmentFile(args.input_bam_file, "rb")

for read in samfile.fetch():
    if not read.is_unmapped and not read.is_secondary and read.is_proper_pair and read.is_paired and\
            not read.is_duplicate and not read.is_supplementary:
        outf.write("@" + read.query_name + '\n')
        outf.write(read.query_sequence + '\n')
        outf.write("+" + '\n')
        outf.write(read.to_string().split('\t')[10] + '\n')

outf.close()
