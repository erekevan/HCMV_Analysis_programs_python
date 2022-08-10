#! /bin/bash

"""
Author: Ahmed Al Qaffas
Program: HCMV Terminal Analyzer
Purpose:
        * To get a, b, and c repeats sequence, frequency, and copy number.
Logic:
A visual representation of the logic is  shown in the following plot link:

'https://drive.google.com/file/d/1bn6IFPEWWAViSuAFGOWZjP4tfyE9Dtwx/view?usp=sharing'

<a> Extract CCS reads that map to a constructed a sequence. # Constructing a sequence.
<ab> Map (a) reads to b sequence. # UL left or right frequency.
    >> Map (ab) to UL left.
    >> Map (ab) to UL right.
<abc> Map (ab) to c sequence. # Junction sequence analysis
<ab_no_c> Remove reads from (ab) that mapped c sequence.
<c> Extract CCS reads that map to c sequence. # Right end c-a c-0 analysis.
    >> Remove all reads that has a_sequence in it.
<c_no_b> Remove all reads that has b_sequence from (c).

<Visualization> Visualization is done by IGV using any of the BAM files

<Analysis> Each repeat has
"""
# Load libraries
import argparse
import os
import sys
import subprocess
import pandas as pd
import numpy as np
from IPython.display import display
from termcolor import colored

# Add minimap2 to path:
sys.path.append("/opt/homebrew/bin/")


## Methods:
# Take every new line as a command to be executed in bash
# In: String
# Out: None
def _run_supressed(s):
    for l in s.split(os.linesep):
        subprocess.run(['bash', '-c', l], capture_output=True, stdout=None, stderr=None)
        # Created by Darren Chan: https://github.com/darrennchan8


'''_________________________________________________________________________________________________________________'''


# Goal: extrac reads within a region mapped the reference
# In: reads file, coords start (s) and end (e), and the name to output
# Out: aligment file to
# def _extract_me(reads='', name='EndStage_', start=1, end=10):
#     p = subprocess.run(['bash', '-c', f"samtools view {reads} | cut -f3 | head -1"],
#                        capture_output=True, stdout=None, stderr=None)
#     ref_name = str(p.stdout, 'utf-8').strip()
#     _run_supressed(f"""
#     samtools view -Sbh {reads} "{ref_name}:1-10" > {name}.bam
#     samtools index {name}.bam
#     """)



'''_________________________________________________________________________________________________________________'''


def _u_complment_me_uno(ref='file_name.fsat', name='some'):
    _run_supressed(f"""
    seqtk seq -r {ref} > {name}_revcom.fa
    """)
    return f'{name}_revcom.fa'


'''_________________________________________________________________________________________________________________'''


# Map reference to reads and index it
# Used for:
# Getting reads that align to only to the reference or without the reference.
# IGV Visualization.
# In: str: ref (fasta), reads (fastq), name, F_value (F 4 or f 4)
#    : int: q_value (quality score)
# Out: str: name (saved as)
# Note: F0x900 allows to keep the best alignment of any given read.
def _mapper(ref, reads, name="AddName_p1", F_value="F 4", q_value=30):
    name = name + ".bam"
    _run_supressed(f"""
    minimap2 -H --secondary=no -ax map-hifi {ref} {reads} | samtools view -F0x900 -Sbh -{F_value} -q{q_value} | samtools sort -o {name} 
    samtools index {name}   
    """)
    print(f'>> Mapping {reads} to {ref}')
    return name


'''_________________________________________________________________________________________________________________'''


# Convert Bam to Fastq
# In: str: bam_file (acquired from _mapper)
# Out: str: fastq file name (saved output)
def _converter(bam_file="AddName_p1"):
    out = bam_file.split('.')[0] + '.fastq'
    print(out, colored("<- filtered reads file name", 'magenta'))
    _run_supressed(f"""
    samtools fastq {bam_file}> {out}
    """)
    print(f">> converting {bam_file} to {out}")
    return out


'''_____________________________________M__ab____________________________________________________________________________'''


# Get the total number of reads mapped to the reference.
# In: str: file name
# Out: str: saved file name
def _total_reads(file_name):
    p = subprocess.run(['bash', '-c', f"samtools view -c -F0x900 {file_name}"],
                       capture_output=True, stdout=None, stderr=None)
    output_dis = str(p.stdout, 'utf-8').strip()
    print(colored(f"The number of reads mapped to '{file_name}' is: ", 'green'), colored(output_dis, 'white'))
    return int(output_dis)


'''_________________________________________________________________________________________________________________'''


# Using a common segment of a file name, running this code:
# Creates a combined  file between two or more fastq files.
# create a text file that has the frequency and header of reads in the combined file.
# In: str: matching_string (target file extension), out_name (name to be saved), header
# Out: int: out_numb (number of reads found
def _mergfastq(matching_string="", out_name="name", header="ccs"):
    _run_supressed(f"cat {matching_string}.fastq > {out_name}.fastq")
    p = subprocess.run(['bash', '-c', f"grep -e \"{header}\" {out_name}.fastq | sort | uniq -c | wc"],
                       capture_output=True)
    subprocess.run(['bash', '-c',
                    f"grep -e \"{header}\" {out_name}.fastq | sort | uniq -c > analysis_outputs/{out_name}_headrs_count.txt"])
    out_numb = str(p.stdout, 'utf-8').strip().split("     ")[0]
    print(f"The merged files using regex {matching_string} and header {header} resulted in\n {out_numb} reads")
    return out_numb


'''_________________________________________________________________________________________________________________'''


# remap all created fastq files to the reference genome.
# In: str: ref (name of the reference file)
#   : list: list of obtained
# Out: None
def _create_vis(ref, los):
    for i in los:
        _mapper(ref, i, name=f'ref__{i.split(".")[0]}')




'''_________________________________________________________________________________________________________________'''


# This code remap using NGML tool (more accurate but computationally more expensive)
# Create VCF files which can be used to determine the variants between the a sequences.
# Normalize the variants.
# In: str: ref, reads
#   : int: threads, depth
# Out: None
# Note: F0x900 keeps best scored alignment of any given read.
def _create_VCF(ref="", reads="", threads=8, depth=8000, name='V__'):
    ref_name = name
    _run_supressed(f"""mkdir analysis_outputs/{ref_name}_analysis""")
    _run_supressed(f"""
    ngmlr -x pacbio -t {threads} -r {ref} -q {reads} | samtools view -Sbh -F0x900 -o ngml_{ref_name}_realn.bam
    samtools sort -o sorted_ngml_{ref_name}_realn.bam ngml_{ref_name}_realn.bam
    samtools index sorted_ngml_{ref_name}_realn.bam
    bcftools mpileup -d {depth} -f {ref} sorted_ngml_{ref_name}_realn.bam | bcftools call -mv -Oz -o ngml_calls_{ref_name}_reds.vcf.gz 
    bcftools index ngml_calls_{ref_name}_reds.vcf.gz
    bcftools norm -f {ref} ngml_calls_{ref_name}_reds.vcf.gz | egrep -v "^#" > out_{ref_name}_variants.txt
    mv *{ref_name}* analysis_outputs/{ref_name}_analysis/
    """)


'''_________________________________________________________________________________________________________________'''


# Create ad output table of frequencies as a webpage (html)
# In: dict: key:: name of repeat, value:: counts
def _create_freq_table(d):
    headers = ['Number of reads', 'Percentage']
    # create matrix
    mat = np.zeros((len(d.keys()), len(headers)))
    outtable = pd.DataFrame(mat, columns=headers, index=d.keys())
    for k, v in d.items():
        v1 = v[0]
        v2 = v[1]
        for i in range(2):
            outtable.loc[k, headers[0]] = v1
            outtable.loc[k, headers[1]] = str(round(v2 * 100, 4)) + " %"
    display(outtable.sort_values(by=['Number of reads']))  # show the table in Terminal
    # generate pretty table
    outtable.to_html('analysis_outputs/table_of_frequencies.html',
                     index=True, header=True, justify='right')


'''_________________________________________________________________________________________________________________'''


def super_fun(ref, name, reads, d, ccs_count, l, exclude=False, F='F 4'):  # Reduce functions
    if exclude:
        mapped = _mapper(ref=ref, name=name, reads=reads, F_value="f 4", q_value=0)
    else:
        mapped = _mapper(ref=ref, name=name, reads=reads, F_value=F)  # Mapping reads to reference
    number_of_reads = _total_reads(mapped)
    d[name] = [number_of_reads, number_of_reads / ccs_count]
    converted = _converter(mapped)
    l.append(converted)
    return mapped, converted


'''_________________________________________________________________________________________________________________'''
# Running the script################ Running the script ############################### Running the script #############
'''_________________________________________________________________________________________________________________'''

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
    # Terminal Region
    parser.add_argument('-ult', '--ult', metavar='ult',
                        help="+2000 bp of UL from the left genomic region in FastA format")
    parser.add_argument('-ulto', '--ulto', metavar='ulto',
                        help="+2000 bp of UL from the left genomic region in FastA format from another isoform")
    # Junction region - UL
    parser.add_argument('-ulj', '--ulj', metavar='ulj',
                        help="+2000 bp of UL from the junction  region in FastA format")
    # Junction isoform
    parser.add_argument('-uljo', '--uljo', metavar='uljo',
                        help="+2000 bp of UL from the junction  region in FastA format")

    # Similarly, for each end of the US terminal region, we need two files that capture each isoform sequence.
    parser.add_argument('-ust', '--ust', metavar='ust',
                        help="+2000 bp of US from the terminal genomic region in FastA format")
    parser.add_argument('-usto', '--usto', metavar='usto',
                        help="+2000 bp of US from the terminal genomic region in FastA format from another isoform")
    # Junction region -US
    parser.add_argument('-usj', '--usj', metavar='usj',
                        help="+2000 bp of US from the junction genomic region in FastA format")
    # Junction isoform
    parser.add_argument('-usjo', '--usjo', metavar='usjo',
                        help="+2000 bp of US from the junction genomic region in FastA format")

    # Parse arguments:
    args = parser.parse_args()
    # initiate some logistics:

    # The following folder will contain the following stats:
    # Number of reads that has been merged from two fastq files and their total count.
    # Stats regrading a, b, and c sequences.
    _run_supressed("mkdir analysis_outputs")

    # list that contains which sequence got generated for visualization.
    los = []

    # dict that contains reads counts. Used to create table of Frequencies.
    freq_dict = {}

    # Starting the analysis pipeline fo a sequence:

    """
    (A) Mapping step:
     each mapping step will generate:
        i- bam alignment file 
        ii-fastq file of filtered reads. 
     * sometimes the mapping step will also be added to a que to be mapped against the whole genomic sequence. 
    """
    if not args.ccs:  # Pipeline requires CCS/PacBio reads
        print("<<<<<<<<<<<You did not use any reads!!>>>>>>>>>>>")
        exit()

    if not args.bseq or not args.cseq or not args.aseq:
        print("Although you can run this program without all the repeats included some of \nthe "
              "analysis steps demands that all the repeats to be presents.")

    ccs_count = 1  # Sanity check
    if args.ccs:
        p1 = subprocess.run(['bash', '-c',
                             f'echo $(cat {args.ccs} | wc -l)/4|bc'], capture_output=True)
        ccs_count = int(str(p1.stdout, 'utf-8').strip())
    print(colored(f"The number of CCS reads is:", 'red'), colored(ccs_count, 'white', attrs=['bold']))
    """
<a>    
    """
    if args.aseq:
        print(colored('<a>', 'red', attrs=['underline', 'bold']))
        # create a prime
        print('creating reverse complement of the a sequence (saved as a_revcom.fa)')
        a_prime = _u_complment_me_uno(ref=args.aseq, name='a')
        # Map CCS reads to a_seq only:
        print('The follwoing reads contains both a and a_inverted')
        a_bam, a_fastq = super_fun(ref=args.aseq, name="M__a", reads=args.ccs, d=freq_dict, l=los,
                                   ccs_count=ccs_count)
        # Create VCF
        _create_VCF(ref=args.aseq, reads=a_fastq, name='V__all_a')

        # ''' a direct'''
        # ad_bam, ad_fastq = super_fun(ref=args.aseq, name="M__ad", reads=args.ccs, d=freq_dict, l=los,
        #                              ccs_count=ccs_count, F='F 16')
        # _create_VCF(ref=args.aseq, reads=ad_fastq, name='V__ad')
        #
        # ''' a reverse'''
        # ar_bam, ar_fastq = super_fun(ref=args.aseq, name="M__ar", reads=args.ccs, d=freq_dict, l=los,
        #                              ccs_count=ccs_count, F='f 16')
        # _create_VCF(ref=a_prime, reads=ar_fastq, name='V__ar')
        '''
<b>    
        '''
    if args.bseq:
        print(colored('<b>', 'red', attrs=['underline', 'bold']))
        # create b prime
        print('creating reverse complement of the b sequence (saved as b_revcom.fa)')
        b_prime = _u_complment_me_uno(ref=args.bseq, name='b')
        # Map CCS reads to b_seq only:
        b_bam, b_fastq = super_fun(ref=args.bseq, name='M__b', reads=args.ccs, d=freq_dict, l=los,
                                   ccs_count=ccs_count)
        # Create VCF
        _create_VCF(ref=args.bseq, reads=b_fastq, name='V__all_b')

        # ''' b direct'''
        # bd_bam, bd_fastq = super_fun(ref=args.bseq, name="M__bd", reads=args.ccs, d=freq_dict, l=los,
        #                              ccs_count=ccs_count, F='F 16')
        # _create_VCF(ref=args.bseq, reads=bd_fastq, name='V__bd')
        # ''' b reverse'''
        # br_bam, br_fastq = super_fun(ref=args.bseq, name="M__br", reads=args.ccs, d=freq_dict, l=los,
        #                              ccs_count=ccs_count, F='f 16')
        # _create_VCF(ref=b_prime, reads=br_fastq, name='V__br')

        if args.aseq:
            print(colored('<ab>', 'red', attrs=['underline', 'bold']))
            """
<ab>
                """
            # Map (a) to b_seq:
            ab_bam, ab_fastq = super_fun(ref=args.bseq, name='M__ab', reads=a_fastq, d=freq_dict, l=los,
                                         ccs_count=ccs_count)
            # # Map (a direct) to b_direct
            # print(colored('<adbd>', 'yellow', attrs=['underline', 'bold']))
            # adbd_bam, adbd_fastq = super_fun(ref=args.bseq, reads=ad_fastq, name='M__adbd', d=freq_dict, l=los,
            #                                  ccs_count=ccs_count)
            # print(colored('<arbr>', 'yellow', attrs=['underline', 'bold']))
            # # Map (a reverse) to b_prime
            # arbr_bam, arbr_fastq = super_fun(ref=b_prime, reads=ar_fastq, name='M__arbr', d=freq_dict, l=los,
            #                                  ccs_count=ccs_count)

        """
<c>
        ! Note: these reads contain both junction and terminal sequences. 
            """
    if args.cseq:
        print(colored('<c>', 'red', attrs=['underline', 'bold']))
        # create c prime
        print('creating reverse complement of the c sequence (saved as c_revcom.fa)')
        c_prime = _u_complment_me_uno(ref=args.cseq, name='c')

        # Map CCS reads to c_seq:
        c_bam, c_fastq = super_fun(ref=args.cseq, name='M__c', reads=args.ccs, d=freq_dict, l=los,
                                   ccs_count=ccs_count)
        _create_VCF(ref=args.cseq, reads=c_fastq, name='V__all_c')

        # ''' c direct'''
        # cd_bam, cd_fastq = super_fun(ref=args.cseq, name="M__cd", reads=args.ccs, d=freq_dict, l=los,
        #                              ccs_count=ccs_count, F='F 16')
        # _create_VCF(ref=args.cseq, reads=c_fastq, name='V__cd')
        #
        # ''' b reverse'''
        # cr_bam, cr_fastq = super_fun(ref=args.cseq, name="M__cr", reads=args.ccs, d=freq_dict, l=los,
        #                              ccs_count=ccs_count, F='f 16')
        # _create_VCF(ref=c_prime, reads=c_fastq, name='V__cr')
        """
<abc>
            """
        if args.aseq and args.bseq:
            print(colored('<abc>', 'red', attrs=['underline', 'bold']))

            abc_bam, abc_fastq = super_fun(ref=args.cseq, name='M__abc', reads=ab_fastq, d=freq_dict,
                                           l=los,
                                           ccs_count=ccs_count)
        """
<a'b'c'>
            """
        # if args.aseq and args.bseq:
        #     print(colored("<a'b'c'>", 'red', attrs=['underline', 'bold']))
        #     abc_prime_bam, abc_prime_fastq = super_fun(ref=c_prime, name='M__abc_prime', reads=arbr_fastq,
        #                                                d=freq_dict, l=los,
        #                                                ccs_count=ccs_count, F='F 16')
        """
<ab_no_c>
            """
        # Filter c_seq from reads containing (ab):
        print(colored('<ab_no_c>', 'red', attrs=['underline', 'bold']))
        ab_no_c_bam, ab_no_c_fastq = super_fun(
            ref=args.cseq, reads=ab_fastq, name="M__ab_no_c", d=freq_dict, l=los, ccs_count=ccs_count,
            exclude=True)

        '''              L Segment Analysis                 '''

        """
<ab_no_c_prime>
            """
        # # Filter c_prime_seq from reads containing (ab):
        # print(colored('<ab_no_c_prime>', 'red', attrs=['underline', 'bold']))
        # ab_no_c_prime_bam, ab_no_c_prime_fastq = super_fun(
        #     ref=c_prime, reads=ab_fastq, name="M__ab_no_c_prime", d=freq_dict, l=los, ccs_count=ccs_count,
        #     exclude=True)

        """
<ab_no_c_UL_P_terminal>    
            """
        # Map reads that has no c sequences to UL left and right
        # UL prototype
        if args.ult:
            print(colored('Creating FRFs for the L segment of the prototype isoform', 'white', attrs=['bold']))
            print(colored('<ab_no_c_UL_P_terminal>', 'red', attrs=['underline', 'bold']))
            ab0c_ULP_bam, ab0c_ULP_fastq = super_fun(
                ref=args.ult, reads=ab_no_c_fastq, name='M__ul_p_ter', d=freq_dict, l=los, ccs_count=ccs_count
            )
            # create variants count of subreads of a & b
            # # a:
            # _create_VCF(ref=args.aseq, reads=ab0c_ULP_fastq, name='V__ab0c_ULP_vs_a_seq')
            # # b:
            # _create_VCF(ref=args.bseq, reads=ab0c_ULP_fastq, name='V__ab0c_ULP_vs_b_seq')
            # Sanity Check 1
            all_R1_bam, all_R1_fastq = super_fun(
                ref=args.ult, reads=args.ccs, name='M__ccs_vs_left_of_L_p', d=freq_dict, l=los,
                ccs_count=ccs_count
            )

            """
<ab_no_c_UL_P_Jun>
                """
            # This set of reads
            # 1- should match UL_IL_Ter (in case of 1:1 prevalence).
            # 2- if hard stop at end of a > IL_ter
            if args.ulj:
                print(colored('<ab_no_c_UL_P_Jun>', 'red', attrs=['underline', 'bold']))

                ab0c_ULP_jun_bam, ab0c_ULP_jun_fastq = super_fun(
                    ref=args.ulj, reads=ab_no_c_fastq, name='M__uljun_p', d=freq_dict, l=los,
                    ccs_count=ccs_count
                )
                all_R150A_bam, all_R150A_fastq = super_fun(
                    ref=args.ulj, reads=args.ccs, name='M__ccs_vs_right_of_L_p', d=freq_dict, l=los,
                    ccs_count=ccs_count
                )
        """
<ab_no_c_UL_IL_terminal>    
            """
        # Map reads that has no c sequences to UL left and right
        # UL prototype
        if args.ulto:
            print(colored('Creating FRFs for the L segment of IL isoform', 'white',
                          attrs=['bold']))
            print(colored('<ab_no_c_UL_IL_terminal>', 'red', attrs=['underline', 'bold']))
            ab0c_ULIL_bam, ab0c_ULIL_fastq = super_fun(
                ref=args.ulto, reads=ab_no_c_fastq, name='M__ul_il_ter', d=freq_dict, l=los,
                ccs_count=ccs_count)
            """
<ab_no_c_UL_IL_Jun>
                """
            if args.uljo:
                print(colored('<ab_no_c_UL_IL_Jun>', 'red', attrs=['underline', 'bold']))

                ab0c_ULIL_jun_bam, ab0c_ULIL_jun_fastq = super_fun(
                    ref=args.uljo, reads=ab_no_c_fastq, name='M__uljun_il', d=freq_dict, l=los,
                    ccs_count=ccs_count
                )
        '''              S Segment Analysis                     '''
        """
<c_no_b>
            """
        print(colored('Creating FRFs for the S segment of P isoform', 'white', attrs=['bold']))
        # Filter b_seq from reads containing (c):
        if args.bseq:
            print(colored('<c_no_b>', 'red', attrs=['underline', 'bold']))

            c_no_b_bam, c_no_b_fastq = super_fun(
                ref=args.bseq, reads=c_fastq, name='M__c_no_b', exclude=True, d=freq_dict, l=los,
                ccs_count=ccs_count
            )
            """
<US_jun_P>        
                """
            if args.usj:
                print(colored('<US_jun_P>', 'red', attrs=['underline', 'bold']))

                US_jun_P_bam, US_jun_P_fastq = super_fun(
                    ref=args.usj, reads=c_no_b_fastq, name='M__usjun_p', d=freq_dict, l=los, ccs_count=ccs_count
                )
                """
<US_jun_IS>
                    """
                if args.usjo:
                    print(colored('<US_jun_IS>', 'red', attrs=['underline', 'bold']))

                    US_jun_IS_bam, US_jun_IS_fastq = super_fun(
                        ref=args.usjo, reads=c_no_b_fastq, name='M__usjun_IS', d=freq_dict, l=los,
                        ccs_count=ccs_count
                    )
                """
<ca_no_b>
                    """
            print(colored('<ca_no_b>', 'red', attrs=['underline', 'bold']))
            ca0b_bam, ca0b_fastq = super_fun(
                ref=args.aseq, reads=c_no_b_fastq, name='M__ca0b', d=freq_dict, l=los, ccs_count=ccs_count
            )

            """
<US_ter_P>        
                """
            if args.ust:
                print(colored('<US_ter_P>', 'red', attrs=['underline', 'bold']))

                US_P_bam, US_P_fastq = super_fun(
                    ref=args.ust, reads=ca0b_fastq, name='M__us_p_ter', d=freq_dict, l=los, ccs_count=ccs_count
                )
                if args.usto:
                    print(colored('<US_ter_IS>', 'red', attrs=['underline', 'bold']))
                    US_IS_bam, US_IS_fastq = super_fun(
                        ref=args.usto, reads=ca0b_fastq, name='M__us_is_ter', d=freq_dict, l=los,
                        ccs_count=ccs_count
                    )

        """
<c_no_a>
            """
        # Filter a_seq from reads containing (c):
        if args.aseq:
            print(colored('<c_no_a>', 'red', attrs=['underline', 'bold']))

            c_no_a_bam, c_no_a_fastq = super_fun(
                ref=args.aseq, reads=c_fastq, name='M__c_no_a', exclude=True, d=freq_dict, l=los,
                ccs_count=ccs_count
            )
            # The following reads are meant to find inner variants within the C terminal sequence

            """
<c_no_a_US_right_P>
                """
            if args.ust:
                print(colored('<c_no_a_US_right_P>', 'red', attrs=['underline', 'bold']))
                # Map reads in P orientation
                c0a_us_ter_P_bam, c0a_us_ter_P_fastq = super_fun(
                    ref=args.ust, reads=c_no_a_fastq, name='M__c0a_us_p_ter', d=freq_dict, l=los,
                    ccs_count=ccs_count
                )

                """
<c_no_a_US_right_IS>
                    """
                if args.usto:
                    print(colored('<c_no_a_US_right_IS>', 'red', attrs=['underline', 'bold']))

                    c0a_us_ter_IS_bam, c0a_us_ter_IS_fastq = super_fun(
                        ref=args.usto, reads=c_no_a_fastq, name='M__c0a_us_is_ter', d=freq_dict, l=los,
                        ccs_count=ccs_count, F='F 16'
                    )

    """
< * vs WGS>
    """
    if args.ref:
        _mapper(ref=args.ref, reads=args.ccs, name="REF_vs_CCS")
        _create_vis(ref=args.ref, los=los)
        _run_supressed("""
        mkdir analysis_outputs/mapped_to_ref/
        mv ref__* REF_vs_* analysis_outputs/mapped_to_ref/ 
        """)
    _create_freq_table(freq_dict)
