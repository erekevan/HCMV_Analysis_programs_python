#! /bin/bash
"""
Author: Ahmed Al Qaffas
Program: Stats Cuts
Summary: the goal of this program is to generate an N-length fasta file from fastq and map it to a reference then
            create a stat file that contains counts of each read starting and ending position. Also, it will create
            a reverse complement fastq file from the reads and re-run the filtering processes.
Inputs:
        1- Reads file in fastq format
        2- Reference file in fasta format
        3- (optional) number of reads to out put from
                a) left end of reads
                b) right end of reads
Outputs:
        1- Sorted alignment file and its index in bam format.
        2- Two Fasta files of the Nth numbers of reads selected
        3- Two text files that include a count of reads start/ending position in with respect to the alignment
        4- Two fasta files of the Nth number of rev-com reads.
        5- Two text files that include a count of rev-comb reads start/ending position in with respect to the alignment
        6- Four plots that show a genomic graph
"""

import argparse
import os
import subprocess


# Take every new line as a command to be executed in bash
def _run_supressed(s):
    for l in s.split(os.linesep):
        subprocess.run(['bash', '-c', l], capture_output=True, stdout=None, stderr=None)


# _converter converts FastQ to FastA as pre-processing step
def _converter(fq_file="AddName_p1"):
    fq_file = fq_file.split('.')[-2]
    _run_supressed(f"""
    seqtk seq -a {fq_file}.fastq > SC__{fq_file}.fasta
    """)
    print(f'SC__{fq_file}.fasta was created')
    return f'SC__{fq_file}.fasta'


# _get_Ns takes FastA file and return the first N nucleotides from each termini
def _get_Ns(reads='some_fasta', left_end=20, right_end=20):
    reads = reads.split('.')[0]
    _run_supressed(f"""
    seqkit subseq -r 1:{left_end} {reads}.fasta >   SC__{reads}_first_{left_end}_bp.fasta
    seqkit subseq -r -{right_end}:-1 {reads}.fasta >  SC__{reads}_last_{right_end}_bp.fasta
    """)
    print(f"getting N's from SC__{reads}_first_{left_end}_bp.fasta",
          f" and SC__{reads}_last_{right_end}_bp.fasta")
    return f"SC__{reads}_first_{left_end}_bp.fasta", f"SC__{reads}_last_{right_end}_bp.fasta"


# _mapper map newly created sub reads with reference
# '-T' minimum score to output [1]
def _mapper(ref='some_ref.fa', reads='some_reads.fa', name=""):
    _run_supressed(f"""
    echo "indexing reference"
    bwa index {ref}
    bwa mem  -T 1 {ref} {reads} | samtools view -Sbh  | samtools sort -o {name}.bam
    samtools index {name}.bam    
    """)

    print(name, ' alignment was created!')


# _get_counts generate a text file that contains two column of either of the following:
#               > frequency of reads terminated at, bp position
##              OR
#               > frequency of reads started at, bp position
def __get_counts(bam_file=""):
    wokring_string = "bedtools bamtobed -i " + \
                     f"{bam_file}" + " | awk -F '\\t' '{print $2}' | sort | uniq -c | sort -nr"
    _run_supressed(f"""
    {wokring_string} > {bam_file.split('.')[0]}.txt
    """)
    # print(f"Running\n{wokring_string}\n")
    print("Getting counts at each position")


# _flip_me_bb takes a fasta file an reverse complement all the reads in it
# it output a new fastq file with rev-combed reads.
def _flip_me_bb(fq=''):
    _run_supressed(f"""
    seqkit seq -r -p {fq} > revcombed_{fq.split('.')[0]}.fastq
    """)
    print(f"Fastq rev-comb file: revcombed_{fq.split('.')[0]}.fastq  was created!")

    return f"revcombed_{fq.split('.')[0]}.fastq"


def _plot_me(file='file_name', output="output", ab=2000, ca=233500, Jstart=194252, Jend=197770):
    _run_supressed(f"""
    HCMV_Analysis_Plots -c {file} -o {output} -ab {ab} -ca {ca} -js {Jstart} -je {Jend} 
    """)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    This program generate multiple files that counts the number of a repeats sequence""")
    # Taking inputs:
    parser.add_argument('-r', '--ref', metavar='ref', help="WGS/ Ref sequence in FastA format")
    parser.add_argument('-fq', '--fastq', metavar='fq', help="Raw reads in fastq format")
    parser.add_argument('-ls', '--left_end', metavar='ls',
                        help=" Reads to be extracted from right [20]")
    parser.add_argument('-rs', '--right_end', metavar='rs',
                        help="Cutoff number on the right end of each read [20]")
    parser.add_argument('-o', '--output', metavar='output', help="output directory/files name")

    # Parse arguments:
    args = parser.parse_args()

    if args.fastq:
        # pre-processing original reads
        reads = _converter(args.fastq)
        # pre-processing rev-comb reads
        revcom = _flip_me_bb(args.fastq)  # rev comb
        rc_reads = _converter(revcom)
        LR, RR = _get_Ns(reads, left_end=args.left_end, right_end=args.right_end)
        rcLR, rcRR = _get_Ns(rc_reads, left_end=args.left_end, right_end=args.right_end)
        if args.ref:
            # Mapping all reads to reference:
            # _mapper(ref=args.ref, reads=reads, name='SC_ref_VS_all_reads')
            # __get_counts('SC_ref_VS_all_reads.bam')
            # _plot_me(file='SC_ref_VS_all_reads.txt', output='ref_VS_all_reads_sc')

            # Getting counts for original reads:
            _mapper(ref=args.ref, reads=LR, name='SC__original_left_end')
            __get_counts('SC__original_left_end.bam')
            _plot_me(file='SC__original_left_end.txt', output='original_left_end')

            _mapper(ref=args.ref, reads=RR, name='SC__original_right_end')
            __get_counts('SC__original_right_end.bam')
            _plot_me(file='SC__original_right_end.txt', output='original_right_end')

            # Getting counts for reverse complement reads:
            _mapper(ref=args.ref, reads=rcLR, name='SC__revcomb_left_end')
            __get_counts('SC__revcomb_left_end.bam')
            _plot_me(file='SC__revcomb_left_end.txt', output='revcomb_left_end')
            _mapper(ref=args.ref, reads=rcRR, name='SC__revcomb_right_end')
            __get_counts('SC__revcomb_right_end.bam')
            _plot_me(file='SC__revcomb_right_end.txt', output='revcomb_right_end')

    # Generate plots
    # Move all files to a new directory
    _run_supressed(f"""
    mkdir stats_cuts_outputs
    mv SC__* stats_cuts_outputs
    """)

    # Run Rscript
