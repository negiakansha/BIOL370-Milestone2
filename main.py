'''
Course: Introduction to Bioinformatics
Assignment: Midterm Milestone #2
Date: Nov 22, 2024

This python3 program is designed to recursively globally align a set
of protein sequences into one alignment.
'''

def ValidateProteins(proteinSeqs):
    print("This function currently does nothing.\n")
    return True

def main():
    #Read input of FASTA formatted protein sequences
    inputFile = open("alignment.txt", "r").readlines()

    #Store the protein sequences along with the names of those sequences
    nameOfSeqs = []
    proteinSeqs = []
    readNextLine = False

    #Sort through input file
    for eachLine in inputFile:
        #Save name for the file
        if readNextLine:
            #Save dna seq associated with that line and find transversions/transitions
            proteinSeqs.append(eachLine.strip())
            readNextLine = False
        if ">" in eachLine:
            #Skip the first char and save name to array if the bounds allow it
            nameOfSeqs.append(eachLine[1 : len(eachLine)].strip())
            #read the next line and store the DNA sequence
            readNextLine = True

    #Validate protein sequences
    validProteinSeq = ValidateProteins(proteinSeqs)
    if not validProteinSeq:
        print("\nA valid protein sequence was not given. Please try again.")
        return
    
    print("Here are the protein sequences given: ")
    print(proteinSeqs)

#run main function
main()