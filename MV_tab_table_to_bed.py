"""
Author: E
Goal: parse MacVector exported feature table and create a BED file.
"""
import sys

# Get reference name
def _getref(file):
    with open(file, 'r') as f:
        mew = f.read()
        mew = mew.split("strain=")
        return mew[1].split("\t")[0].split('\n')[0]


try:
    if len(sys.argv) > 1:
        fileIn = sys.argv[1]
        print(fileIn)
    else:
        fileIn = input("Please enter features file name:")

    refName = ""
    if not refName:
        refName = input("We couldn't locate reference name, please input type it in: ")
    with open(fileIn, "r") as f:
        to_out = {}
        with open('out.bed', 'w') as outfile:
            for line in f:
                start, end, flipped, name = "", "", "", ""
                if line.startswith("CDS"):
                    print(line)
                    line = line.split("\t")
                    start = line[1].strip()
                    end = line[2].strip()
                    flipped = "+"
                    if line[3] == "C":
                        flipped = "-"
                    name = line[4].split("product=")[1].split("/")[0].strip()
                elif line.startswith("repeat_region"):
                    line = line.split("\t")
                    start = line[1].strip()
                    end = line[2].strip()
                    name = line[4].strip().split(' ')[0].split("=")[1]
                if len(start) > 0:
                    if ' ' in start:
                        print(f"gene {name} has more than one exon")
                        for i in range(len(start.split(' '))):
                            s = start.split(' ')[i]
                            e = end.split(' ')[i]
                            outfile.write("{}\t{}\t{}\t{}\t{}\t{}".format(refName, s, e, name, 0, flipped))
                            outfile.write('\n')
                    else:
                        to_out[name] = [refName, start, end, name, 0, flipped]
                        outfile.write("{}\t{}\t{}\t{}\t{}\t{}".format(refName, start, end, name, 0, flipped))
                        outfile.write('\n')
        print("The output file saved under the name 'out.bed'")
except FileNotFoundError:
    print(f"File {fileIn} does not exist... now exiting")
    exit(0)

