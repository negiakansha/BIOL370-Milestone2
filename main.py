'''
Course: Introduction to Bioinformatics
Assignment: Midterm Milestone #2
Date: Nov 22, 2024

This python3 program is designed to recursively globally align a set
of protein sequences into one alignment.
'''

#Make sure there are at least two protein sequences in protein sequence input
def ValidateProteins(proteinSeqs):
    if len(proteinSeqs) > 1:
        return True
    return False

def GloballyAlign(protein1, protein2):
    # ????
    return True

def AlignSequences(proteinSeqs):
    # Analyze a set of protein sequences to find which pair shows highest similarity
    # Return the two sequences with the most matching protein sequences
    similarity = 0
    tempSimilarity = 0
    highestMatchIndex1 = -1
    highestMatchIndex2 = -1

    for i in range(len(proteinSeqs)):
        for j in range(i + 1, len(proteinSeqs)):
            protein1 = proteinSeqs[i]
            protein2 = proteinSeqs[j]
            if len(protein1) < len(protein2):
                smallerProteinLen = len(protein1)
            else:
                smallerProteinLen = len(protein2)

            similarity = 0  # Reset similarity
            for charIndex in range(smallerProteinLen):
                if protein1[charIndex] == protein2[charIndex]:
                    similarity += 1

            if similarity > tempSimilarity:
                highestMatchIndex1 = i
                highestMatchIndex2 = j
                tempSimilarity = similarity  # Store highest similarity so far in tempSimilarity

    print([highestMatchIndex1, highestMatchIndex2])
    print("num matches: " + str(tempSimilarity))

    # globally align those two sequences using your team’s substitution matrix
    #GloballyAlign(highestMatchIndex1, highestMatchIndex2)

    # Replace the two protein sequences with their alignment in the original set of sequences
    newProteinSeqs = []
    alignmentPlaced = False
    for i in range(len(proteinSeqs)):
        if not alignmentPlaced and (i == highestMatchIndex1 or i == highestMatchIndex2):
            newProteinSeqs.append("replace here") # GloballyAlign(highestMatchIndex1, highestMatchIndex2)
            alignmentPlaced = True
        elif i != highestMatchIndex1 and i != highestMatchIndex2:
            newProteinSeqs.append(proteinSeqs[i])

    print(newProteinSeqs)
    print(len(newProteinSeqs))
    proteinSeqs = newProteinSeqs

    return newProteinSeqs

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
    print(len(proteinSeqs))

    #re-iterate the steps until all sequences have been included in a single multiple alignment. 
    # while len(proteinSeqs) > 1:
    proteinSeqs = AlignSequences(proteinSeqs)
    AlignSequences(proteinSeqs)

#run main function
main()