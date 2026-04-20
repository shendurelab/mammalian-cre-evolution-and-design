# # Opens sequencing files then outputs list of barocdes and ORFs
# # @author Diego 07/22/2021
# # modified by JB 04/03/2022: 
# added all in/out as paths as arguments
# pull additional information from bam file (alignment flag, positions-->serves as "UMI")
# changed to fixed sequence position search after barcode

from Bio import SeqIO
import sys, gzip, pysam, re, argparse


parser = argparse.ArgumentParser('.')
parser.add_argument('--input_bam', '-i', help='Position sorted BAM (or list of bams), from bowtie to list of regions.')
parser.add_argument('--input_BC1_fastq', '-ifq1', help='fastq.gz file corresponding to oBC read.')
parser.add_argument('--output_file', '-o', help='Tab delimited file with cell, mutation barcode, read count, umi count. All observed barcodes correctable to a whitelist are reported.')
parser.add_argument('--barcode1_start', type=int, help='start position of barcode in read')
parser.add_argument('--barcode1_end', type=int, help='length of barcode')
parser.add_argument('--search_seq1', help='Sequence to search for immediately after barcode.')

args = parser.parse_args()

input_bam=args.input_bam
input_BC1_fastq=args.input_BC1_fastq

output_file=args.output_file

barcode1_end=args.barcode1_end
barcode1_start=args.barcode1_start
search_seq1=args.search_seq1
search_seq1_start = barcode1_start + barcode1_end + 1 
search_seq1_end = barcode1_start + barcode1_end + 1 + len(search_seq1)

print('reading in ' + input_BC1_fastq)

# to deal with multi-mapper, easier to START by reading in the BC fastq and create a READID to BC dictionary. 
read_BC1=dict()
with gzip.open(input_BC1_fastq, 'rt') as BC1_fq:
    bh = SeqIO.parse(BC1_fq, 'fastq') 
    for read in bh:
        seq = read.seq.upper()
                
        #print(seq)
        #print(seq[search_seq_start: search_seq_end])
        
        # searching for the constant sequence after the BC for valid read
        is_ok=False
        if search_seq1=="no_seq":
            is_ok=True
        elif seq[search_seq1_start: search_seq1_end]==search_seq1:
            is_ok=True
        
        if is_ok:
            barcode1 = seq[barcode1_start: barcode1_start+barcode1_end+1]
            read_BC1[read.id]=barcode1
            #print(barcode)
            #print(read_alignments[read.id])
            #out.write("\t".join([str(read.id), str(barcode), str(read_alignments[read.id])]) + '\n')
print("done reading in BC1 fastq file")
  



print('reading in ' + input_bam)

# extracting some information from the bam file
with gzip.open(output_file, 'wt') as out:

    out.write("\t".join(['read_id', 'BC1', 'CRE_id', 'mapq', 'CRE_read', 'read_forward_in_ref','is_read1', 'read_start', 'read_end', 'align_length']) + '\n')
    
    #read_alignments=dict()
    with pysam.AlignmentFile(input_bam, "rb") as samfile:
        for read in samfile:
            # note, each paired-end read will be present twice in the output file, read1 and read2
            # restrict to looping to mate pairs with is_read1?
            if read.reference_name is not None and read.mapping_quality is not None and read.query_name in read_BC1:                    
                summary="\t".join([read.reference_name, str(read.mapping_quality), str(read.query_sequence), str(read.is_reverse), str(read.is_read1), str(read.reference_start), str(read.reference_end), str(read.reference_length)])
                #read_alignments[read.query_name]=summary
                out.write("\t".join([str(read.query_name), str(read_BC1[read.query_name]), str(summary)]) + '\n')
                    











