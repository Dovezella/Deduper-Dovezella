# Author: Dove Enicks

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.5"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = set("ACTGNactgn")
RNA_bases = set("AUGCNaugcn")

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    phred_score = ord(letter)
    qual_score = phred_score -33
    return qual_score

def qual_score(phred_score: str) -> float:
    """this will calculate an average quality score from the whole string"""
    sum = 0
    for value in phred_score:
        q_score=convert_phred(value)
        sum+= q_score
    return sum/len(phred_score)

def validate_base_seq(sequence: str,RNAflag=False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs, and Ns. False otherwise. Case insensitive.'''
    seq = set(sequence.upper())
    return seq<=(RNA_bases if RNAflag else DNA_bases) 

def gc_content(sequence: str) -> float:
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(sequence), "String contains invalid characters - are you sure you used a DNA sequence?"
    
    DNA = sequence.upper()
    Gs = DNA.count("G")
    Cs = DNA.count("C")
    return (Gs+Cs)/len(DNA)

def calc_median(lst: list) -> float:
    '''This will calculate the median value at each position in the sequence for a file'''
    if len(lst)%2 == 0:
        c=lst[(len(lst))//2-1]
        b=lst[(len(lst))//2]
        a=(c+b)/2
    else:
        a=(lst[(len(lst))//2])
                
    return a

def oneline_fasta(filein: str, fileout: str) -> str:
    '''This will convert a multiline fasta file to a single line file (sequences are all on one line wihtout new line characters). Don't forget to specify filein and fileout when you call it'''
    with open (filein, "r") as fin, open(fileout,"w") as fout:
        line = fin.readline()
        fout.write(line)
        for line in fin: 
                if line [0] == ">":
                    fout.write(f'\n{line}')
                else:
                    fout.write(line.strip())
    return

rev_dict={"A":"T","T":"A","G":"C","C":"G","N":"N"}

def rev_comp(seq: str, rev_dict) -> str:
    '''This function will take a sequence and return the reverse complement of that sequence'''
    assert bioinfo.validate_base_seq(seq)
    seq=seq[::-1]
    rev=""
    for i in seq:
        rev += rev_dict[i]
    return rev

def read_rec(fq_filehand) -> list:          #This function is important to be able to read files R1-R4 one record at a time, which I used later in a while true loop
    '''This function can read the first record of the file called with the file handle. It should be used in a while true loop to be able to go through all the records and can be used in conjuntion with other files at the same time by calling them by their file handles'''
    header=fq_filehand.readline().strip()
    seq=fq_filehand.readline().strip()
    plus=fq_filehand.readline().strip()
    quality=fq_filehand.readline().strip()
    return [header,seq,plus,quality]

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")
    #tests for qual_score
    assert qual_score("EEE") == 36
    assert qual_score("#I") == 21
    assert qual_score("EJ") == 38.5
    print("You calcluated the correct average phred score")
    #tests for validate_base_seq
    assert validate_base_seq("AATAGAT", False) == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("TATUC",False) == False
    assert validate_base_seq("UCUGCU", False) == False
    print("Passed DNA and RNA tests")
    #tests for gc_content
    assert gc_content("GCGCGC") == 1
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATGCAT") == 0.5
    print("correctly calculated GC content")
    #need to make assert tests for calc_median
    assert calc_median([1,2,3,4,5]) == 3
    assert calc_median([1,2,3,4,5,6]) == 3.5
    
