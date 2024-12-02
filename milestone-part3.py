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

#This is the substitution matrix we are using to determine alignment for protein sequences
BLOSUM62 = {
    'A': {'A': 4, 'C': 0, 'D': -2, 'E': -1, 'F': -2, 'G': 0, 'H': -1, 'I': -1, 'K': -1, 'L': -1, 'M': -1, 'N': -2, 'P': -1, 'Q': -1, 'R': -1, 'S': 1, 'T': 0, 'V': 0, 'W': -3, 'Y': -2},
    'C': {'A': 0, 'C': 9, 'D': -3, 'E': -4, 'F': -2, 'G': -3, 'H': -3, 'I': -1, 'K': -3, 'L': -1, 'M': -1, 'N': -3, 'P': -2, 'Q': -3, 'R': -3, 'S': -1, 'T': -1, 'V': -1, 'W': -2, 'Y': -2},
    'D': {'A': -2, 'C': -3, 'D': 6, 'E': 2, 'F': -3, 'G': -1, 'H': 0, 'I': -3, 'K': -1, 'L': -3, 'M': -2, 'N': 1, 'P': -1, 'Q': 0, 'R': 0, 'S': 0, 'T': -1, 'V': -2, 'W': -4, 'Y': -3},
    'E': {'A': -1, 'C': -4, 'D': 2, 'E': 5, 'F': -3, 'G': -2, 'H': 0, 'I': -3, 'K': 1, 'L': -2, 'M': -2, 'N': 0, 'P': -1, 'Q': 2, 'R': 1, 'S': 0, 'T': -1, 'V': -2, 'W': -3, 'Y': -2},
    'F': {'A': -2, 'C': -2, 'D': -3, 'E': -3, 'F': 6, 'G': -3, 'H': -1, 'I': 0, 'K': -2, 'L': 0, 'M': 0, 'N': -3, 'P': -2, 'Q': -2, 'R': -2, 'S': -1, 'T': -1, 'V': 0, 'W': 1, 'Y': 3},
    'G': {'A': 0, 'C': -3, 'D': -1, 'E': -2, 'F': -3, 'G': 6, 'H': -2, 'I': -4, 'K': -2, 'L': -4, 'M': -3, 'N': -1, 'P': -2, 'Q': -2, 'R': -2, 'S': 0, 'T': -2, 'V': -3, 'W': -2, 'Y': -3},
    'H': {'A': -1, 'C': -3, 'D': 0, 'E': 0, 'F': -1, 'G': -2, 'H': 8, 'I': -2, 'K': -1, 'L': -2, 'M': -2, 'N': 1, 'P': -1, 'Q': 0, 'R': 0, 'S': -1, 'T': -2, 'V': -2, 'W': -1, 'Y': 2},
    'I': {'A': -1, 'C': -1, 'D': -3, 'E': -3, 'F': 0, 'G': -4, 'H': -2, 'I': 4, 'K': -1, 'L': 2, 'M': 1, 'N': -3, 'P': -2, 'Q': -1, 'R': -1, 'S': -2, 'T': -1, 'V': 3, 'W': -3, 'Y': -1},
    'K': {'A': -1, 'C': -3, 'D': -1, 'E': 1, 'F': -2, 'G': -2, 'H': -1, 'I': -1, 'K': 5, 'L': -1, 'M': -1, 'N': 0, 'P': -1, 'Q': 1, 'R': 2, 'S': 0, 'T': -1, 'V': -2, 'W': -3, 'Y': -2},
    'L': {'A': -1, 'C': -1, 'D': -3, 'E': -2, 'F': 0, 'G': -4, 'H': -2, 'I': 2, 'K': -1, 'L': 4, 'M': 2, 'N': -3, 'P': -2, 'Q': -1, 'R': -1, 'S': -2, 'T': -1, 'V': 1, 'W': -2, 'Y': -1},
    'M': {'A': -1, 'C': -1, 'D': -2, 'E': -2, 'F': 0, 'G': -3, 'H': -2, 'I': 1, 'K': -1, 'L': 2, 'M': 5, 'N': -2, 'P': -2, 'Q': 0, 'R': 0, 'S': -1, 'T': -1, 'V': 1, 'W': -1, 'Y': -1},
    'N': {'A': -2, 'C': -3, 'D': 1, 'E': 0, 'F': -3, 'G': -1, 'H': 1, 'I': -3, 'K': 0, 'L': -3, 'M': -2, 'N': 6, 'P': -1, 'Q': 0, 'R': 0, 'S': 0, 'T': -1, 'V': -3, 'W': -4, 'Y': -3},
    'P': {'A': -1, 'C': -2, 'D': -1, 'E': -1, 'F': -2, 'G': -2, 'H': -1, 'I': -2, 'K': -1, 'L': -2, 'M': -2, 'N': -1, 'P': 7, 'Q': -1, 'R': -1, 'S': -1, 'T': -1, 'V': -2, 'W': -3, 'Y': -2},
    'Q': {'A': -1, 'C': -3, 'D': 0, 'E': 2, 'F': -2, 'G': -2, 'H': 0, 'I': -1, 'K': 1, 'L': -1, 'M': 0, 'N': 0, 'P': -1, 'Q': 5, 'R': 1, 'S': 0, 'T': -1, 'V': -2, 'W': -3, 'Y': -2},
    'R': {'A': -1, 'C': -3, 'D': 0, 'E': 1, 'F': -2, 'G': -2, 'H': 0, 'I': -1, 'K': 2, 'L': -1, 'M': 0, 'N': 0, 'P': -1, 'Q': 1, 'R': 5, 'S': 0, 'T': -1, 'V': -2, 'W': -3, 'Y': -2},
    'S': {'A': 1, 'C': -1, 'D': 0, 'E': 0, 'F': -1, 'G': 0, 'H': -1, 'I': -2, 'K': 0, 'L': -2, 'M': -1, 'N': 0, 'P': -1, 'Q': 0, 'R': 0, 'S': 4, 'T': 1, 'V': -2, 'W': -3, 'Y': -2},
    'T': {'A': 0, 'C': -1, 'D': -1, 'E': -1, 'F': -1, 'G': -2, 'H': -2, 'I': -1, 'K': -1, 'L': -1, 'M': -1, 'N': -1, 'P': -1, 'Q': -1, 'R': -1, 'S': 1, 'T': 5, 'V': 1, 'W': -2, 'Y': -2},
    'V': {'A': 0, 'C': -1, 'D': -2, 'E': -2, 'F': 0, 'G': -3, 'H': -2, 'I': 3, 'K': -2, 'L': 1, 'M': 1, 'N': -3, 'P': -2, 'Q': -2, 'R': -2, 'S': -2, 'T': 1, 'V': 4, 'W': -3, 'Y': -1},
    'W': {'A': -3, 'C': -2, 'D': -4, 'E': -3, 'F': 1, 'G': -2, 'H': -1, 'I': -3, 'K': -3, 'L': -2, 'M': -1, 'N': -4, 'P': -3, 'Q': -3, 'R': -3, 'S': -3, 'T': -2, 'V': -3, 'W': 11, 'Y': 2},
    'Y': {'A': -2, 'C': -2, 'D': -3, 'E': -2, 'F': 3, 'G': -3, 'H': 2, 'I': -1, 'K': -2, 'L': -1, 'M': -1, 'N': -3, 'P': -2, 'Q': -2, 'R': -2, 'S': -2, 'T': -2, 'V': -1, 'W': 2, 'Y': 7}
}

gap_penalty = -2

# def GloballyAlign(protein1, protein2):
#     # ????
#     return True

def GloballyAlign(protein1, protein2):
    #Initialize matrix based on protein1 length and protein2 length
    matrixWidth = len(protein1)
    matrixHeight = len(protein2)
    matrix = [[0 for x in range(matrixHeight+1)] for y in range(matrixWidth+1)]
    
    #Edit cells in first row and column based on gap penalty
    for i in range(matrixWidth+1):
        matrix[i][0] = i * gap_penalty
    for j in range(matrixHeight+1):
        matrix[0][j] = j * gap_penalty
    
    # Edit the rest of the matrix cells with max value calculated
    for i in range(1, matrixWidth+1):
        for j in range(1, matrixHeight+1):            
            # For each cell (i, j), calculate:
            # (a) Cell (i – 1, j) – Gap penalty
            # (b) Cell (i, j – 1) – Gap penalty
            # (c) Cell (i – 1, j – 1) ± Match/mismatch
            a = matrix[i - 1][j] + gap_penalty
            b = matrix[i][j - 1] + gap_penalty
            # find match or mismatch
            matchValue = BLOSUM62[protein1[i-1]].get(protein2[j-1])
            c = (matrix[i - 1][j - 1]) + matchValue

            # Then fill in the largest of these values
            matrix[i][j] = max(a, b, c)
    
    # Perform the backward pass (trace the optimal alignment)
    aligned_protein1 = []
    aligned_protein2 = []
    
    i = matrixWidth
    j = matrixHeight
    while i > 0 and j > 0:
        score = matrix[i][j]
        score_diag = matrix[i-1][j-1]
        score_up = matrix[i-1][j]
        score_left = matrix[i][j-1]
        
        # If the score came from the diagonal, it's a match/mismatch
        if score == score_diag + BLOSUM62[protein1[i-1]].get(protein2[j-1]):
            aligned_protein1.append(protein1[i-1])
            aligned_protein2.append(protein2[j-1])
            i -= 1
            j -= 1
        # If the score came from the left, it's an insertion in protein1
        elif score == score_left + gap_penalty:
            aligned_protein1.append('-')
            aligned_protein2.append(protein2[j-1])
            j -= 1
        # If the score came from the up, it's a deletion in protein2
        else:
            aligned_protein1.append(protein1[i-1])
            aligned_protein2.append('-')
            i -= 1
    
    # If we've reached the first column or row, fill with gaps
    while i > 0:
        aligned_protein1.append(protein1[i-1])
        aligned_protein2.append('-')
        i -= 1
    
    while j > 0:
        aligned_protein1.append('-')
        aligned_protein2.append(protein2[j-1])
        j -= 1
    
    # Reverse the alignments to get them in correct order
    aligned_protein1 = ''.join(reversed(aligned_protein1))
    aligned_protein2 = ''.join(reversed(aligned_protein2))
    
    return aligned_protein1, aligned_protein2


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

    print("matches are at indices: ")
    print([highestMatchIndex1, highestMatchIndex2])
    print("num matches: " + str(tempSimilarity))

    # globally align those two sequences using your team’s substitution matrix
    #GloballyAlign(highestMatchIndex1, highestMatchIndex2)

    # Replace the two protein sequences with their alignment in the original set of sequences
    newProteinSeqs = []
    alignmentPlaced = False
    for i in range(len(proteinSeqs)):
        if not alignmentPlaced and (i == highestMatchIndex1 or i == highestMatchIndex2):
            #REPLACE HERE WITH GLOBALLY ALIGNED OR ALIGNED STRING
            newProteinSeqs.append(proteinSeqs[i]) # GloballyAlign(highestMatchIndex1, highestMatchIndex2)

            # Testing the function with two example protein sequences
            protein1 = "PLEASANT"
            protein2 = "MEANINGFUL"

            aligned1, aligned2 = GloballyAlign(protein1, protein2)
            print(f"Protein 1 aligned: {aligned1}")
            print(f"Protein 2 aligned: {aligned2}")

            alignmentPlaced = True
        elif i != highestMatchIndex1 and i != highestMatchIndex2:
            newProteinSeqs.append(proteinSeqs[i])

    print(newProteinSeqs)
    print("Length of sequences: " + str(len(newProteinSeqs)))
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
    print("Length of sequences: " + str(len(proteinSeqs)))

    #re-iterate the steps until all sequences have been included in a single multiple alignment. 
    while len(proteinSeqs) > 1:
        proteinSeqs = AlignSequences(proteinSeqs)

#run main function
main()