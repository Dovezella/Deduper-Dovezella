# Deduper: Pseudo-code

## The problem:  
Due to amplification during sequencing, there can be over-representation of different reads in a sample due to various errors in sequencing or sometimes some sequences get amplified more than others. Since some sequences are inherently not as common as others in the original sample, after amplification these differences might be covered up or not as obvious especially if we just want to be looking at one unique read instead of multiple copies of the same read, which especially important for expression data in transcriptomics. So for this program, we are trying to write code which is able to parse through a sam file to find unique sequences based off their chromosome, corrected starting position, umi, and strandedness to write out a sam file which contains only one of each of those unique reads.

#### Please reference my test files here:
```/home/dovee/bioinfo/Bi624/Deduper-Dovezella/unit_test```

### Pseudo-Functions needed:
1.
rev_comp function from my bioinfo.py which will reverse complement a given sequence

2.
def UMI_list(input known UMI file: string) -> output known_UMI[list of known UMI's]
    ```This function will parse through a file that contains known UMI and return a list of them```

expected input: STL96.txt
expected output: [list of 96 known UMI's from given file]

3.
**(note:  are the input/output correct to say file: string or how would you approach that?)**
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

4.
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

<!-- **note: DONT DO THIS, but it is important to consider how to loop through for the deduping and also consider errors in umi**
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
call pos_correct on sam file
use samtools on corrected_sam to create a sorted sam file based off UMI, start pos, chromosome, and strand
then runs python script on corrected, sorted file to output deduped sam file
delete intermediate file

#### dedupe python pseudo code:

open new file to write 
open position corrected, sorted sam file
read.line on sam file
write headers to new file (using @ symbol)
<!-- initialize f_position_count{key = position, value = number of reads with that position } 
initialize r_position_count{key = position, value = number of reads with that position } -->
when at first read split line by tabs and remove new line character
call bitwise_strand function on line
save UMI, chromosome, start pos, and bitwise_strand(T or F) in potential_list and whole line as potential_line
go to next line and begin for line loop:
    strip, and tab separate, and call bitwise_strand function on line -> store in compare_list and store line in compare_line
    if potential_list /= compare_list in any of the four components:
        if forward strand:
            change corrected start pos back to original by add appropriate nt's from soft clip or no change if not clipped
            write out potential_line to new sam file
            make compare_list and compare_line the new potential_line and potential_list
        elif reverse
            revert back to original start pos and UMI
            write out potential_line to new sam file
            make compare_list and compare_line the new potential_line and potential_list
    elif potential_list == compare_list:
        continue to next line and rewrite compare_list and compare_line with next line info
    elif at end of file
        close files

## Consider:

Making a report to output number of reads excluded and percentages of kept reads.
Function for fixing start position  




     


