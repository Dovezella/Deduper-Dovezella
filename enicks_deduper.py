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
    '''This will parse the bitwise flag of a read in SAM file for its bitwise flag and decode it to determine strandedness.
    Strip and split a readline and pass this function the 3rd position of the list for the bitwise flag'''
    if ((int(flag) & 16) == 16): 
        rev = "rev"
    else:
        rev = "forward"

    return rev     

def adjust_5_position (line: list) -> str:
    '''This function will take in the current readline (stripped and split) of a SAM file, take out the reported 5-prime, 
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
    mathm = 0
    if bit == "forward":
        Lclip = CIGAR[0:4]
        if "S" in Lclip:  #only left-handed soft clipping
            math = re.search("^[0-9]+", Lclip)
            math = int(math.group(0))          
            f_adjusted = int(unadjusted) - math 
    elif bit == "rev":                                                
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
        if "M" in CIGAR:
            mathm = re.findall("[0-9]+(?=M)", CIGAR)
            mathm = list(map(int, mathm))
            mathm = sum(mathm)
    r_adjusted = rmath + dmath + nmath + mathm + int(unadjusted) -1
    
    return str(f_adjusted) if bit == "forward" else str(r_adjusted)

        

out = open(outfile, "w")    #open file to write unique reads 
dupe = open(dupe, "w")      #open file to write duplicates
file = open(filename, "r")  #open input sorted, sam file
err = open(erunk, "w")      #open file to write unknown umi reads or error umi reads

my_UMIs = Umi_list(umi)     #populate set of all known reads

compare_reads={}    #dictionary for comparing reads where {key = list[UMI, strand, adjusted 5' position], value = count of duplicates}
header=""           #initialize header variable, so you can call line after all headers have been read and not skip first line
for header in file:         
    if header.startswith("@"):             #write out the headers
        out.write(header)
    else:
        break           

line_2_write = header #grab first read line from last line taken from header loop which should be first read line
line = line_2_write.strip().split()
current_chrom = line[2]                 #assign chromosome for line comparisons, to know when to start over with empty dictionary
split_header = line[0].split(":")
current_read_umi = split_header[7]      #assign umi of current read
if current_read_umi in my_UMIs:
    out.write(line_2_write)             #write out first read line to unique as it is the first so not a duplicate
    read_strand = bitwise_strand(line[1])       #assign strandedness
    start_pos = adjust_5_position(line)         #adujst start position
    compare_list= current_read_umi + read_strand + start_pos    #create string of comparison values for dictionary
    compare_reads[compare_list]= 0                              #assign values for line comparison of first read line to dictionary
else:
    err.write(line_2_write)         #else, this is not a known umi

while True:                             #Begin loop for comparisons to find duplicates
    compare_2_write = file.readline()   #next readline is saved to write 
    compare = compare_2_write.strip().split()   
    if compare_2_write == "":                   #should only break if end of file
        break
    elif compare[2] != current_chrom:           #if new chromosome, begin over with line comparisons 
        current_chrom = compare[2]     #assign new chromosome for line comparisons
        compare_reads={}                        #empty dictionary as duplicates can not share chromosome
        out.write(compare_2_write)              #write out last line of previous chromosome as not a duplicate
        compare_split_header_umi = compare[0].split(":")
        compare_read_umi = compare_split_header_umi[7]  #assign umi
        strand = bitwise_strand(compare[1])             #assign strandedness
        five_start = adjust_5_position(compare)         #assign adjusted position
        comparison = compare_read_umi + strand + five_start     #set variable to compare these values
        compare_reads[comparison] = 0                   #add variable to dictionary as unique
        line_2_write = file.readline()                  #this is where you grab the second line of new chromosome 
        line = line_2_write.strip().split()
        # current_chrom = line[2]     #assign new chromosome for line comparisons
        if current_chrom == line[2]:              #if chromosome is the same, we need to compare to previous
            split_header = line[0].split(":")
            current_read_umi = split_header[7]      #assign umi for next line
            if current_read_umi in my_UMIs:
                strand = bitwise_strand(line[1])
                five_start = adjust_5_position(line)
                comparison = current_read_umi + strand + five_start
                if comparison in compare_reads:
                    dupe.write(line_2_write)
                    compare_reads[comparison]+=1
                else:
                    out.write(line_2_write)
                    compare_reads[comparison]=0
            else:
                err.write(line_2_write)
        else: # current_chrom != line[2]: #this line is new chromosome first line, so add to dictionary three values = 0 and move on   #so we want to write it out
            current_chrom = line[2]
            compare_reads={}
            split_header = line[0].split(":")
            current_read_umi = split_header[7]      #assign umi for next line
            if current_read_umi in my_UMIs:         #check if umi is known
                out.write(line_2_write)             #write out as first line should not be duplicate
                read_strand = bitwise_strand(line[1])   #assign strandedness
                start_pos = adjust_5_position(line)     #assign adjusted position
                compare_list=current_read_umi + read_strand + start_pos #combine as values to add to comparison dictionary
                compare_reads[compare_list]= 0          #assign variable to dictionary
            else:
                err.write(line_2_write)                 #if not known, write to error-unknown file
    else:       #new line chromosome = current_chrom
        compare_split_header_umi = compare[0].split(":")
        compare_read_umi = compare_split_header_umi[7]  #assign umi
        if compare_read_umi in my_UMIs:                 #check if umi is known
            strand = bitwise_strand(compare[1])         #assign strandedness
            five_start = adjust_5_position(compare)     #assign adjusted position
            comparison = compare_read_umi + strand + five_start     #combine values as variable to add to comparison dictionary
            if comparison in compare_reads:             #check if variable already in comparions dictionary
                dupe.write(compare_2_write)             #if so, write out to duplicate file
                compare_reads[comparison]+= 1           #count up number of duplicates
            else:
                out.write(compare_2_write)              #if not already in dictionary, write out line to unique sam file
                compare_reads[comparison] = 0           #add to comparison dictionary to find duplicates
        else:
            err.write(compare_2_write)                  #if not known umi, write to error-unknown file
    


out.close()
dupe.close()
file.close()
err.close()         #close files :)