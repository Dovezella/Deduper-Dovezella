# Deduper: Pseudo-code

## The problem:  
Due to amplification during sequencing, there can be over-representation of different reads in a sample due to various errors in sequencing or sometimes some sequences get amplified more than others. Since some sequences are inherently not as common as others in the original sample, after amplification these differences might be covered up or not as obvious especially if we just want to be looking at one unique read instead of multiple copies of the same read, which especially important for expression data in transcriptomics. So for this program, we are trying to write code which is able to parse through a sam file to find unique sequences based off their chromosome, corrected starting position, umi, and strandedness to write out a sam file which contains only one of each of those unique reads.

#### Please reference my test files here:
```/home/dovee/bioinfo/Bi624/Deduper-Dovezella/unit_test```

### Pseudo-Functions needed:
1.
def UMI_list(input known UMI file: string) -> output known_UMI[list of known UMI's]
    ```This function will parse through a file that contains known UMI and return a list of them```

expected input: STL96.txt
expected output: [list of 96 known UMI's from given file]

2.
def position_correct(input ordered sam file: string) -> output new, corrected sam file:
    ```This function will take an ordered same file, and go line by line to correct the starting position of each read based off of the bit-wise flag and using the CIGAR string. It will also reverse complement the UMI of reverse strands and exclude UMI's with errors in sequence or unknown UMI's```
    open UMI file to get known UMI's and call UMI_list on it
    open new file to write to
    read.line of ordered sam file
    write header lines to output sam file using @'s
    get read
    for each line [check column 2 (bitwise flag) for strandedness and if UMI is known]:
        line strip and split by tab
        if N in UMI, 
            do not write out line and continue to next line
        elif UMI not in UMI_list,
            do not write out line and continue to next line
        else call bitwise_strand on line 
            if forward, 
                check cigar string (column 6) for soft clipping
                if clipped at beginning, 
                    subtract number of clipped nt to start position
                    write out to new sam file
            elif reverse,
                rev_comp the UMI 
                check cigar string
                if no clipping at beginning,
                    add total of nt's minus 1 in cigar string to start pos  
                    write out to new sam file  
                if clipped at beginning,
                    subtract clip number +1 from read length (71) and add to start pos (or just add numbers after initial soft clip and add to start pos)
                    write out to new sam file
    close files
    
    return corrected_sam

expected input:  /home/dovee/bioinfo/Bi624/Deduper-Dovezella/unit_test/pos_correct_test_input.sam
expected output:  /home/dovee/bioinfo/Bi624/Deduper-Dovezella/unit_test/unit_test_output.sam

3.
def bitwise_strand(input list position of bitwise flag column in sam file) -> output True for Forward, False for Reverse:
    ```This will take a line that has been split by tabs from a sam file, and given the bitwise flag column it will parse for the strandedness. It will return a boolean for whether it is forward (T) or reverse (F) which can be evaluated downstream```

expected input:
    Forward:
        16 bit is not set (0)
    Reverese:
        16 bit is set (1) 
expected output:
    Forward: 
        True
    Reverse:
        False


## Then to actually complete my assignment:

#### Samtools on initial sam file (be sure to sort)
samtools view -b file.sam > file.bam
samtools sort file.bam -o file.sorted.bam
samtools view -h file.sorted.bam > file.sorted.sam

#### dedupe python pseudo code:

Pass the sam file a sorted by chromosome and start position 

open UMI file
open duplicate file to write to
open desired deduped output sam file to write to

initialize set for known umi's
read file to populate set of known umi's

initialize dicitonary for chromosome-based line comparisons {key = list[UMI, strand, adjusted 5' position], value = count of duplicates}

Write out the headers to the output deduper file

Read the first line,
split line by tab and strip new line character
assign chromosome to var "current_chrom"
if UMI is known
    write out to deduped sam file
    decode bitwise flag for strand 
    decode CIGAR string to adujst 5' position
    add all three elements to dictionary and count set to 0
else
    unknown, write out to duplicate file and continue
while true 
Read second line,
If chromosome of new line /= current_chrom,
    empty dictionary
    Read the first line,
    assign chromosome to var "current_chrom"
    if UMI is known
        write out to deduped sam file
        decode bitwise flag for strand 
        decode CIGAR string to adujst 5' position
        add all three elements to dictionary and count set to 0
    else
        unknown, write out to duplicate file and continue
else (new line chromosome = current_chrom)
    If UMI is known 
        decode bitwise for strand
        decode CIGAR for adjusted 5' position
        if these three values are already in line comparison dictionary
            write out to duplicate file, up counter, and continue to next line
        else three values are not in dictionary,
            add to dictionary and count set to 0 
    else
        unknown, write out to duplicate file and continue

Continue with this process until you get to the end of the file

if readline = empty string/end of file
    break
