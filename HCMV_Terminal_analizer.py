#! /bin/bash

"""
Author: Ahmed Al Qaffas
Program: HCMV Terminal Analyzer
Purpose:
        * To get a, b, and c repeats sequence, frequency, and copy number.
Logic:

<a> Extract CCS reads that map to a constructed a sequence. # Constructing a sequence.
<ab> Map (a) reads to b sequence. # UL left or right frequency.
    >> Map (ab) to UL left.
    >> Map (ab) to UL right.
<abc> Map (ab) to c sequence. # Junction sequence analysis
<ab_no_c> Remove reads from (ab) that mapped c sequence.
<c> Extract CCS reads that map to c sequence. # Right end c-a c-0 analysis.
    >> Remove all reads that has a_sequence in it.
<c_no_b> Remove all reads that has b_sequence from (c).

<Visualization>

"""

import argparse
import os
import subprocess


## Methods:
# Take every new line as a command to be executed in bash
def _run_supressed(s):
    for l in s.split(os.linesep):
        subprocess.run(['bash', '-c', l], capture_output=True, stdout=None, stderr=None)


# Map reference to reads and index it
def _mapper(ref, reads, name="AddName_p1", F_value="F 4", q_value=30):
    _run_supressed(f"""
    minimap2 -H -ax map-pb {ref} {reads} | samtools view -Sbh -{F_value} -q{q_value} | samtools sort -o {name}.bam
    samtools index {name}.bam    
    """)


# Convert Bam to Fastq
def _converter(bam_file="AddName_p1"):
    _run_supressed(f"""
    samtools fastq {bam_file}.bam > {bam_file}.fastq
    """)
    return f'{bam_file}.fastq'


def _total_reads(file_name):
    p = subprocess.run(['bash', '-c', f"samtools view -c {file_name}.bam"],
                       capture_output=True, stdout=None, stderr=None)
    output_dis = f"The number of reads mapped to '{file_name}' is: " + str(p.stdout, 'utf-8').strip()
    print(output_dis)
    return output_dis


def _mergfastq(matching_string="", out_name="name", header="ccs"):
    _run_supressed(f"cat {matching_string}.fastq > {out_name}.fastq")
    p = subprocess.run(['bash', '-c', f"grep -e \"{header}\" {out_name}.fastq | sort | uniq -c | wc"],
                       capture_output=True)
    subprocess.run(['bash', '-c',
                    f"grep -e \"{header}\" {out_name}.fastq | sort | uniq -c > analysis_outputs/{out_name}_headrs_count.txt"])
    out_numb = str(p.stdout, 'utf-8').strip().split("     ")[0]
    print(f"The merged files using regex {matching_string} and header {header} resulted in\n {out_numb} reads")
    return out_numb


def _create_vis(ref, los):
    for i in los:
        _mapper(ref, i, name=f'ref_{i.split(".")[0]}')


# Running the program:
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    This program generate multiple files to estimate the terminal repeats sequence""")
    # Taking inputs:
    parser.add_argument('-a', '--aseq', metavar='aseq', help="Complete 'a' sequence in FastA format")
    parser.add_argument('-b', '--bseq', metavar='bseq', help="Complete 'b' sequence in FastA format")
    parser.add_argument('-c', '--cseq', metavar='cseq', help="Complete 'c' sequence in FastA format")
    parser.add_argument('-ccs', '--ccs', metavar='ccs', help="CCS reads in FastQ format")
    parser.add_argument('-ref', '--ref', metavar='ref', help="Whole genome reference in FastA")

    # For each end of the UL region, we need two files that capture each isoform sequence.
    parser.add_argument('-ult', '--ult', metavar='ult',
                        help="+2000 bp of UL from the left genomic region in FastA format")
    parser.add_argument('-ulto', '--ulto', metavar='ulto',
                        help="+2000 bp of UL from the left genomic region in FastA format from another isoform")

    parser.add_argument('-ulj', '--ulj', metavar='ulj', help="+2000 bp of UL from the junction  region in FastA format")

    # Similarly, for each end of the US region, we need two files that capture each isoform sequence.
    parser.add_argument('-usr', '--ust', metavar='ust',
                        help="+2000 bp of US from the terminal genomic region in FastA format")
    parser.add_argument('-usj', '--usj', metavar='usj',
                        help="+2000 bp of US from the terminal genomic region in FastA format")
    # Parse arguments:
    args = parser.parse_args()

    # initiate some logistics:

    # The following folder will contain the following stats:
    # Number of reads that has been merged from two fastq files and their total count.
    # Stats regrading a, b, and c sequences.
    _run_supressed("mkdir analysis_outputs")

    # list that contains which sequence got generated for visualization.
    los = []

    # Starting the analysis pipeline fo a sequence:
    """
    (A) Mapping step:
     each mapping step will generate:
        i- bam alignment file 
        ii-fastq file of filtered reads. 
     * sometimes the mapping step will also be added to a que to be mapped against the whole genomic sequence. 
    """
    if not args.bseq or not args.cseq or not args.aseq:
        print("Although you can run this program without all the repeats included some of \nthe "
              "steps demands that all the repeats to be presents.")

    if args.aseq:
        """<a>"""
        # Map CCS reads to a_seq only:
        _mapper(ref=args.aseq, name="a", reads=args.ccs)
        a_mapped = _total_reads("a")
        los.append(_converter("a"))

    if args.bseq:
        """<b>"""
        # Map CCS reads to a_seq only:
        _mapper(ref=args.bseq, name="b", reads=args.ccs)
        b_mapped = _total_reads("b")
        los.append(_converter("b"))

        if args.bseq:
            """<ab>"""
            # Map (a) to b_seq:
            _mapper(ref=args.bseq, reads="a.fastq", name="ab")
            ab_mapped = _total_reads("ab")
            los.append(_converter("ab"))

    if args.cseq:
        """<c>"""
        # Map CCS reads to c_seq:
        _mapper(ref=args.cseq, reads=args.ccs, name="c")
        c_mapped = _total_reads("c")
        los.append(_converter("c"))

        if args.aseq and args.bseq:
            """<abc>"""
            # Map (ab) to c_seq:
            _mapper(ref=args.cseq, reads="ab.fastq", name="abc")
            abc_mapped = _total_reads("abc")
            los.append(_converter("abc"))

            """<ab_no_c>"""
            # Filter c_seq from reads containing (ab):
            _mapper(ref=args.cseq, reads="ab.fastq", name="ab_no_c", F_value="f 4", q_value=0)
            ab0c_mapped = _total_reads("ab_no_c")
            los.append(_converter("ab_no_c"))

            # Map reads that has no c sequences to UL left and right
            if args.ult:
                _mapper(ref=args.ult, reads='ab_no_c.fastq', name='ab0c_ul')
                ab0c_UL = _total_reads("ab0c_ul")
                _converter("ab0c_ul")
                if args.ulto:
                    _mapper(ref=args.ulto, reads='ab_no_c.fastq', name='ab0c_ul_iso')
                    ab0c_UL = _total_reads("ab0c_ul_iso")
                    _converter("ab0c_ul_iso")
                    _mergfastq(matching_string="*_ul*", out_name="UL_two_iso")

        if args.aseq:
            """<c_no_a>"""  # remove
            # Filter a_seq from reads containing (c):
            _mapper(ref=args.aseq, reads="c.fastq", name="c_no_a", F_value="f 4", q_value=0)
            print(_total_reads("c_no_a"))
            los.append(_converter("c_no_a"))
        if args.bseq:
            """c_no_b"""
            _mapper(ref=args.bseq, reads="c.fastq", name="c_no_b", F_value="f 4", q_value=0)
            print(_total_reads("c_no_b"))
            los.append(_converter("c_no_b"))
        if args.usj:
            """c-US-junction"""
            _mapper(ref=args.usj, reads="c.fastq", name="c_US_jun")
            print(_total_reads("c_US_jun"))
            los.append(_converter("c_US_jun"))
        if args.ust:
            """US_c terminal"""
            _mapper(ref=args.ust, reads="c.fastq", name="c_US_ter")
            print(_total_reads("c_US_ter"))
            los.append(_converter("c_US_ter"))
    if args.ref:
        print(f"The entered reference  is {args.ref}")
        _create_vis(args.ref, los)

    print(los)
    print("Some calculations and outputs saved on analysis")

