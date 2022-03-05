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
        for line in f.readlines():
            id = ''
            if line.find("CDS"):
                line = line.split('\t')

        print("done")



except FileNotFoundError:
    print(f"File {fileIn} does not exist... now exiting")
    exit(0)
