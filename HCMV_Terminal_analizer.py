#!  /bin/bash
import argparse
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
# To import:
import sys
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


def _total_reads(file_name):
    p = subprocess.run(['bash', '-c', f"samtools view -c {file_name}.bam"],
                       capture_output=True, stdout=None, stderr=None)
    return f"The number of reads mapped to '{file_name}' is: " + str(p.stdout, 'utf-8')


def _create_vis(los):
    for i in ['a']:
        _run_supressed("""
    """)


if __name__ == '__main__':
    print("")
    print("Running analysis:\n")
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', help="Complete 'a' sequence")
    parser.add_argument('-b', help="Complete 'b' sequence")
    parser.add_argument('-c', help="Complete 'c' sequence")
    parser.add_argument('-ccs', help="CCS reads")
    parser.add_argument('-ref', help="Whole genome reference in FastA")
    parser.add_argument('-ult', help="+2000 bp of UL from the left genomic region")
    parser.add_argument('-ulj', help="+2000 bp of UL from the junction  region")
    parser.add_argument('-usr', help="+2000 bp of US from the right genomic region")


    # ref = sys.argv[1]
    # name = ref.split("/")[-1].split(".")[0]
    ccs_reads = sys.argv[1]
    a_seq = sys.argv[2]  # Complete 'a' sequence
    b_seq = sys.argv[3]  # Complete 'b' sequence
    c_seq = sys.argv[4]  # Complete 'c' or c' sequence
    UL_left = sys.argv[5]  # About 3000 bases of  UL near the terminal region
    UL_right = sys.argv[6]  # About 3000 bases of  UL junction region


    """<a>"""
    # Map CCS reads to a_seq only:
    _mapper(ref=a_seq, name="a", reads=ccs_reads)
    print(_total_reads("a"))
    _converter("a")

    """<ab>"""
    # Map (a) to b_seq:
    _mapper(ref=b_seq, reads="a.fastq", name="ab")
    print(_total_reads("ab"))
    _converter("ab")

    """<ab-UL-Left>"""
    # Map (ab) to UL_left:
    _mapper(ref=UL_left, reads="ab.fastq", name="ab_UL_left")
    print(_total_reads("ab_UL_left"))
    _converter("ab_UL_left")

    """<ab-UL-right>"""
    # Map (ab) to UL_left:
    _mapper(ref=UL_right, reads="ab.fastq", name="ab_UL_right")
    print(_total_reads("ab_UL_right"))
    _converter("ab_UL_right")

    """<abc>"""
    # Map (ab) to c_seq:
    _mapper(ref=c_seq, reads="ab.fastq", name="abc")
    print(_total_reads("abc"))
    _converter("abc")

    """<ab_no_c>"""
    # Filter c_seq from reads containing (ab):
    _mapper(ref=c_seq, reads="ab.fastq", name="ab_no_c", F_value="f 4", q_value=0)
    print(_total_reads("ab_no_c"))
    _converter("ab_no_c")

    """<c>"""
    # Map CCS reads to c_seq:
    _mapper(ref=c_seq, reads=ccs_reads, name="c")
    print(_total_reads("c"))
    _converter("c")

    """<c_no_a>"""  # remove
    # Filter a_seq from reads containing (c):
    _mapper(ref=a_seq, reads="c.fastq", name="c_no_a", F_value="f 4", q_value=0)
    print(_total_reads("c_no_a"))
    _converter("c_no_a")

    """c_no_a_no_b"""
    # Filter a_seq from reads containing (c):
    _mapper(ref=b_seq, reads="c.fastq", name="c_no_b", F_value="f 4", q_value=0)
    print(_total_reads("c_no_b"))
    _converter("c_no_b")

    """ Creating visualizations """

    """ Generating Stats"""


