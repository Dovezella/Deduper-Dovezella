#!/usr/bin/env python

import argparse
import re

def get_args():
    parser = argparse.ArgumentParser(description="A program for PS6")
    parser.add_argument("-f", "--filename", help="filename", required=True, type=str)
    parser.add_argument("-o", "--outfile", help="output file name", required= True, type=str)
    parser.add_argument("-u", "--umi", help="designates file containing the list of known UMIs", required= True, type=str)
    parser.add_argument("-d", "--dupes", help="file to create which holds non-unique reads aka duplicates", type = str)
    parser.add_argument("-e", "--extra", help="file to create which holds reads that have erros in UMI or are not known UMIs", type=str)
    return parser.parse_args()

args=get_args()
filename=args.filename
outfile=args.outfile
umi=args.umi
dupe=args.dupes
erunk=args.extra

def Umi_list (filein: str) -> set[str]:
    '''This will read a file that contains only lines, where each line is a unique, known UMI for our deduping process'''
    UMIs=set()
    with open (filein, "r") as fin:
        for line in fin:
            line=line.strip()
            UMIs.add(line)
    
    return UMIs

def bitwise_strand (flag: str) -> str:
    '''This will parse the stripped and split line of a SAM file for its bitwise flag and decode it to determine strandedness.
    True is Reverse, False is Forward'''
    if ((int(flag) & 16) == 16): 
        rev = "+"
    else:
        rev = "-"

    return rev     

def adjust_5_position (line: list) -> str:
    '''This function will take in the current readline of a SAM file, take out the reported 5-prime, 
    left most mapping position, parse the CIGAR string, and output a variable that is the adjusted 
    position based off any relevant CIGAR information'''
    unadjusted = line[3]
    CIGAR = line[5]
    bit = bitwise_strand(line[1])
    f_adjusted = 0 
    r_adjusted = 0
    rmath = 0
    dmath = 0
    nmath = 0
    if bit == "+":
        Lclip = CIGAR[0:4]
        if "S" in Lclip:  #only left-handed soft clipping
            math = re.search("^[0-9]+", Lclip)
            math = int(math.group(0))           #idk what this group error is?
            f_adjusted = int(unadjusted) - math 
    elif bit == "-":                                                
        if CIGAR.endswith("S"):    #only right-handed soft clipping
            rmath = re.search("[0-9]+(?=S$)", CIGAR)
            rmath = int(rmath.group(0))
        if "D" in CIGAR:
            dmath = re.findall("[0-9]+(?=D)", CIGAR)
            dmath = list(map(int, dmath))
            dmath = sum(dmath)
        if "N" in CIGAR:
            nmath = re.findall("[0-9]+(?=N)", CIGAR)
            nmath = list(map(int, nmath))
            nmath = sum(nmath)
    r_adjusted = rmath + dmath + nmath + int(unadjusted) -1
    
    return str(f_adjusted) if bit == "+" else str(r_adjusted)

        

out = open(outfile, "w")
dupe = open(dupe, "w")
file = open(filename, "r")
err = open(erunk, "w")

my_UMIs = Umi_list(umi)

compare_reads={}    #dictionary for comparing reads where {key = list[UMI, strand, adjusted 5' position], value = count of duplicates}
header=""
for header in file:         
    if header.startswith("@"):             #write out the headers
        out.write(header)
    else:
        break           

line_2_write = header #grab first read line from last line taken from header loop which should be first read line
line = line_2_write.strip().split()
current_chrom = line[2]                 #assign chromosome for line comparisons, to know when to start over with empty dictionary
split_header = line[0].split(":")
current_read_umi = split_header[7]
if current_read_umi in my_UMIs:
    out.write(line_2_write)    
    read_strand = bitwise_strand(line[1])
    start_pos = adjust_5_position(line)
    compare_list= current_read_umi + read_strand + start_pos
    compare_reads[compare_list]= 0
else:
    err.write(line_2_write)
# counter = 1     #just for testing
while True:                             #I think the very first read line is getting skipped
    # counter +=1     #just to test
    compare_2_write = file.readline()
    compare = compare_2_write.strip().split()
    if compare_2_write == "": 
        # print("oopsie", f'{counter=}')         #testing, so this shows to be breaking early
        break
    elif compare[2] != current_chrom:
        compare_reads={}
        out.write(compare_2_write)
        compare_split_header_umi = compare[0].split(":")
        compare_read_umi = compare_split_header_umi[7]
        strand = bitwise_strand(compare[1])
        five_start = adjust_5_position(compare)
        comparison = compare_read_umi + strand + five_start
        compare_reads[comparison] = 0
        line_2_write = file.readline()
        line = line_2_write.strip().split()
        current_chrom = line[2]     #assign chromosome for line comparisons, to know when to start over with empty dictionary
        split_header = line[0].split(":")
        current_read_umi = split_header[7]
        if current_read_umi in my_UMIs:
            # print("line 122, new first line", f'{counter=}')    #testing
            out.write(line_2_write)    
            read_strand = bitwise_strand(line[1])
            start_pos = adjust_5_position(line)
            compare_list=current_read_umi + read_strand + start_pos
            compare_reads[compare_list]= 0
        else:
            # print("line 129, unknown", f'{counter=}')   #testing
            err.write(line_2_write)    
    else:       #new line chromosome = current_chrom
        compare_split_header_umi = compare[0].split(":")
        compare_read_umi = compare_split_header_umi[7]
        if compare_read_umi in my_UMIs:
            strand = bitwise_strand(compare[1])
            five_start = adjust_5_position(compare)
            comparison = compare_read_umi + strand + five_start
            if comparison in compare_reads:
                # print("line 129, dupe", f'{counter=}')   #testing
                dupe.write(compare_2_write)
                compare_reads[comparison]+= 1
            else:
                # print("line142, unique", f'{counter=}') #testing
                out.write(compare_2_write)
                compare_reads[comparison] = 0
        else:
            # print("line 149, unknown", f'{counter=}')   #testing
            err.write(compare_2_write)
    


out.close()
dupe.close()
file.close()
err.close()