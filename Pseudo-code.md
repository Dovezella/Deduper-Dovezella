# Deduper: Pseudo-code

## The problem:  
Due to amplification during sequencing, there can be over-representation of different reads in a sample due to various errors in sequencing or sometimes some sequences get amplified more than others. Since some sequences are inherently not as common as others in the original sample, after amplification these differences might be covered up or not as obvious especially if we just want to be looking at one unique read instead of multiple copies of the same read, which especially important for expression data in transcriptomics. So for this program, we are trying to write code which is able to parse through a sam file to find unique sequences based off their chromosome, corrected starting position, umi, and strandedness to write out a sam file which contains only one of each of those unique reads.

#### Please reference my test files here:
```/home/dovee/bioinfo/Bi624/Deduper-Dovezella/unit_test```

### Pseudo-Functions needed:
1.
**(note:  are the input/output correct to say file: string or how would you approach that?)**
def position_correct(input ordered sam file: string) -> output new, corrected sam file:
    ```This function will take an ordered same file, and go line by line to correct the starting position of each read based off of the bit-wise flag and using the CIGAR string```
    open new file to write to
    read.line of ordered sam file
    skip header lines
    get read
    line by line check column 2 (bitwise flag) for strandedness 
    if forward, 
        check cigar string (column 6) for soft clipping
        if clipped at beginning, 
            add number of clipped nt to start position
            write out to new sam file
    if reverse,
        rev_comp the UMI 
        check cigar string
        if there is an insertion,
            add to 71 and add to start pos
        if there is an deletion,
            subtract from 71 and add to start pos    
        if clipped at end,
            subtract clip number from read length (71) and add to start pos
        if skipped portion,
            add to 71 and add to start pos
    write out to new sam file
    close files
    
    return corrected_sam

expected input:  /home/dovee/bioinfo/Bi624/Deduper-Dovezella/unit_test/pos_correct_test_input.sam
expected output:  /home/dovee/bioinfo/Bi624/Deduper-Dovezella/unit_test/unit_test_output.sam

<!-- 2.**note: DONT DO THIS, but it is important to consider how to loop through for the deduping and also consider errors in umi**
def loop_criteria(input ordered sam file: string) -> dictionary{key = umi, value = chromosome}:
    ```This function will take an ordered sam file, and go line by line of the reads to create a dictionary of all the umi's for each chromosome```
    init loop{dictionary}
    open and read sam file line by line
    skip headers
    get first read
    if N in umi,
        skip
    elif umi not in loop{},
        add umi with chromosome as value
    elif umi in loop{},
        if chromosome value not in loop.values{},
            
        if chromosome value in loop.values{},
        continue
    close files
    return loop{}

expected input: 
expected output:  -->


## Then to actually complete my assignment:

#### Slurm Script which:
use samtools to create a sorted sam file 
then runs python script to dedupe the sorted same file

#### dedupe python pseudo code:

call position_correct on sorted sam file
open new file to write unique reads
read.line on position corrected sam file
remove if secondary alignment or unmapped
go back to top of file and read back at top read line
remove if unknown or error in UMI
go back to top of file and read 
write headers to new file
when at first read
initialize f_position_count{key = position, value = number of reads with that position } 
initialize r_position_count{key = position, value = number of reads with that position }
for each umi and chromosome pair,
    check bitwise flag for strandedness 
    if forward:   
        if position start not in f_position_count{},
            add start bp as key and count =1
            write to new output file
        elif position start in f_position_count{},
            continue
    elif line doesn't match umi and chromosome,
        continue  
    if reverse:    
        if position start not in r_position_count{},
            add start bp as key and count =1
            write to new output file
        elif position start in r_position_count{},
            continue
    elif line doesn't match umi and chromosome,
        continue 
close files






     


