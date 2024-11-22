"""
This program gives us log-odd caluclations. 

This praogram takes a text file as an input of protein sequence arranged in FASTA format . It performs several analysis on amino acid sequence. 
There are major 4 functions(along with few helper function): 

1. Amino Acid frequency calculation: 
This function calculates the total frequency of amino acid and individual frequent of amino acid and calculates the proportion. 

2. Observed Freqyenct of Amino Acid Paidrs: 
This function computes the frequency of pairs of amino acid that appear at the sampe position across all sequences. 

3. Expected Frequencies: 
This function calculates frequencies for each pair of amino acids are calculated based on the hypothesis of independent occurrence. This is done by multiplying the individual probabilities (frequencies) of the amino acids.

4. Log odds calculation: 
Using the observed and expected frequencies, the program computes the log-odds ratio for each pair of amino acids. 

5. Output:
    - The program outputs the results in two formats:
        1. Text File (log_odds_matrix.txt)**: Contains the log-odds values for each pair of amino acids in a readable format.
        2. CSV File (log_odds_matrix.csv)**: A matrix of log-odds values that allows easy visualization and further analysis of amino acid pair frequencies.


How to Use:
1. Provide a FASTA format protein sequence file eg: sequence.txt as input. You can change the name on main function.
2. Run the program to calculate the frequencies, observed pair frequencies, expected frequencies, and log-odds scores.
3. Output files will be generated, containing the results in both text and CSV formats.

"""

import math
import csv

#Creating a list for valid amino acid including '-'
valid_amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']


#Function to open the FASTA format txt file to get the sequence
def open_file(file_name):
    sequence_names = []
    sequences = []
    current_sequence = ""   
    try:
        with open(file_name, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith('>'):  # Sequence header line
                    if current_sequence:  # Saving the current sequence before starting a new one
                        sequences.append(current_sequence)
                        current_sequence = ""
                    sequence_names.append(line[1:])  # Saving the sequence name, excluding the '>' character
                else:
                    current_sequence += line  # Append line to the current sequence
            if current_sequence: 
                sequences.append(current_sequence)
    except IOError:
        print("Error: file does not appear to exist")
    return sequence_names, sequences

#Function to validate if there are any invalid characters in the sequences
def validate_aminoacid(sequences, valid_amino_acids):
    sequence_counter = 1  # Initialize the sequence counter, so program could point out the number of sequence which is invalid
    for seq in sequences:  # Iterate through each sequence
        if not all(char in valid_amino_acids for char in seq.upper()):  # Check the individual sequence
            print(f"Invalid amino acid sequence found in sequence number {sequence_counter}")
            exit()  # Terminate the program
        sequence_counter += 1  # Increment the counter for the next sequence

#QNO. 2a 
#This function calculates the frequency of single amino acid and total number of amino acid
def calculate_frequency(sequences):
    #Initialize dictionary for easy counting of individual counting of amino acid 
    counts_for_frequency = {
    'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0,
    'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0,
    'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0, '-': 0 
    }

    total_count = 0

    #Runs for the number of sequence present in the txt file 
    for seq in sequences:
        for char in seq:  # Check if the character is valid
            if char in counts_for_frequency:
                counts_for_frequency[char] += 1  # Increment the count
                total_count += 1  # Increment the total count
    
    # Calculate frequencies
    frequencies = {}
    # Loops through all amino acids in the dictionary
    for amino_acid, count in counts_for_frequency.items():  
        frequencies[amino_acid] = count / total_count

    #Returns to main function 
    return counts_for_frequency, frequencies

#QNO.2b
#This function is set up to calculate the number of amino acid pair 
def calculate_amino_pair(sequences):
    
    # Initialize dictionary for counting occurrences of each amino acid
    counts_for_pair = {
        'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0,
        'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0,
        'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0, '-': 0 
        }

    # List of possible amino acids
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-']

    # Get the length of sequences and the length of individual sequences
    length_of_list = len(sequences)
    length_of_individual_seq = len(sequences[0])

    # Iterate over each position (0 to n-1, where n is the length of the sequences)
    for i in range(length_of_individual_seq):  # i is the position in the sequence
        # Reset the counts for the current position (i)
        for aa in counts_for_pair:
            counts_for_pair[aa] = 0
        
        # Count the occurrences of amino acids at position i
        for sequence in sequences:
            # Access the character at position i (for each sequence)
            if sequence[i] in counts_for_pair:
                counts_for_pair[sequence[i]] += 1  # Increment the count
        
        #Calling the helper function that calculates the pairwise value
        final_pooling = calculate_pairwise(counts_for_pair)  
    
    #This function helps to calculate the proportion
    observed_frequency = calculate_proportion(final_pooling)

    return observed_frequency
    
#Helper function for 2b 
#This function calculates pairs of amino aicd
def calculate_pairwise(counts_for_pair):
    # List of possible amino acids
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-']
    
    #Initializing empty dictionary to store the pairs 
    amino_with_count = {}
    pairwise_frequencies = {}
    
    # Initialize pairwise_frequencies dictionary with all possible pairs set to 0
    for i in range(len(amino_acids)):
        for j in range(i, len(amino_acids)):  # Only counting (i, j) once (no reverse pairs like ('C', 'A'))
            pairwise_frequencies[(amino_acids[i], amino_acids[j])] = 0 + 0.5

    # Populating amino_with_count dictionary with non-zero counts
    for value in counts_for_pair:
        if counts_for_pair[value] != 0:
            amino_with_count[value] = counts_for_pair[value]  # Add the amino acid and its count

    # Loop through the dictionary and calculate pairwise frequencies
    for i, amino_acid1 in enumerate(amino_with_count):  # First loop (outer loop)
        count_i = amino_with_count[amino_acid1]
        
        # Now loop from `i` to the end to form pairs (i, j)
        for j, amino_acid2 in enumerate(list(amino_with_count.keys())[i:]):  # Second loop (inner loop)
            count_j = amino_with_count[amino_acid2]
            
            # If same amino acid, calculate using the self-pair formula
            if amino_acid1 == amino_acid2:  
                pair_count = (count_i * (count_i - 1)) / 2
            else:  # If different amino acids, calculate using the formula: count_i * count_j
                pair_count = count_i * count_j
            
            # Accumulate the frequency into the pairwise_frequencies dictionary
            pairwise_frequencies[(amino_acid1, amino_acid2)] += pair_count

    return pairwise_frequencies

#Helper function for 2b 
#This function calculates the proportion for observed frequency of amino acid pairs 
def calculate_proportion(final_pooling):
    #Calculate the total number of pairs
    total_pairs = sum(final_pooling.values())  # Sum all values in the final_pooling dictionary
    
    # Creating a new dictionary to store the proportions
    proportion_pooling = {}
    
    # Step 3: Calculate the proportion for each pair
    for pair, observed_freq in final_pooling.items():
        proportion = observed_freq / total_pairs if total_pairs > 0 else 0
        proportion_pooling[pair] = proportion

    return proportion_pooling


    #for 
#Qno. 2c
# This function helps to calculate the expected frequency for each amino acid pairs 
def expected_frequencies(frequencies, observed_frequency_of_pairs):
    # Initializing the dictionary to store frequencies
    expected_frequencies_dict = {}

    # Iterate through the observed frequency of pairs
    for (aminoacid1, aminoacid2), observed_freq in observed_frequency_of_pairs.items(): 
        # Get the probability of the first and second amino acids
        Probability_of_aminoacid1 = frequencies[aminoacid1]
        Probability_of_aminoacid2 = frequencies[aminoacid2]

        # Checking if the amino acids are the same or different and calculating accordingly
        if aminoacid1 == aminoacid2:
            # If amino acids are the same
            expected_freq = Probability_of_aminoacid1 * Probability_of_aminoacid1 
        else:
            # If amino acids are different
            expected_freq = Probability_of_aminoacid1 * Probability_of_aminoacid2 * 2
        
        # Store the calculated expected frequency
        expected_frequencies_dict[(aminoacid1, aminoacid2)] = expected_freq

    return expected_frequencies_dict


#Part 2d 
#This function calculates the logodd scores of amino acid pairs 
def Logodd_scores(expected_frequencies_dict, observed_frequency_of_pairs):
    log_odds_dict = {}

    # Loop through the observed frequency pairs to calculate the log odds ratio for each pair
    for (aminoacid1, aminoacid2), observed_freq in observed_frequency_of_pairs.items():
        # Get the expected frequency for this pair
        expected_freq = expected_frequencies_dict.get((aminoacid1, aminoacid2), 0)

        # Avoiding division by zero by checking if the expected frequency is zero, assigning default value of 0 
        if expected_freq == 0:
            log_odds_dict[(aminoacid1, aminoacid2)] = 0
        else:
            # Calculating the log odds ratio
            log_odds = math.log2(observed_freq / expected_freq)
            log_odds_dict[(aminoacid1, aminoacid2)] = log_odds

    return log_odds_dict    

#This function outputs the final matrix in csv and text file 
def log_odds_to_text_and_csv(log_odds_dict, filename_txt='log_odds_matrix.txt', filename_csv='log_odds_matrix.csv'):
    # Write to the text file
    with open(filename_txt, 'w') as file:
        for key, value in log_odds_dict.items():
            # Write each key-value pair to the file, formatted as "key: value"
            file.write(f"{key}: {value}\n")

    # Write to the CSV file
    # Manually create the list of unique amino acids from the dictionary keys
    amino_acids = sorted(set([aa for pair in log_odds_dict.keys() for aa in pair]))
    
    # Open the CSV file in write mode
    with open(filename_csv, mode='w', newline='') as file:
        writer = csv.writer(file)

        #Write the header row (amino acids as column labels)
        writer.writerow([''] + amino_acids)  # First cell is made empty to put in row labels

        #Write each row for each amino acid
        for aa1 in amino_acids:
            row = [aa1]  # Start the row with the amino acid label
            for aa2 in amino_acids:
                # If the pair exists in the dictionary, add the value; otherwise, add 0
                if (aa1, aa2) in log_odds_dict:
                    row.append(log_odds_dict[(aa1, aa2)])
                else:
                    row.append(0)  # Add 0 if the pair does not exist
            writer.writerow(row)

"""
Main function to run the program logic.
"""
def main():
    file_name = "sequence.txt"  #File to read

    #Function to open the files and it passes sequence with name
    #name is not used, but stored incase if needed for futher use
    sequence_names, sequences = open_file(file_name)

    #calling the function that validates aminoacid sequences, teminates the program if invalid amino acid sequence
    validate_aminoacid(sequences, valid_amino_acids)

    #Calulates the indidual and total amino acid frequency 
    counts, individual_amino_frequencies = calculate_frequency(sequences)

    #Calculates the amio acid pair from the sequences 
    observed_frequency_of_pairs = calculate_amino_pair(sequences)

    #Calculates the expected frequency based on individual and total amino acid 
    expected_frequencies_dict = expected_frequencies(individual_amino_frequencies, observed_frequency_of_pairs)

    #Final function that caclulates the log score for the final matrix 
    log_score = Logodd_scores(expected_frequencies_dict, observed_frequency_of_pairs)
    
    #Name for output files, (Change the output filename as needed)
    filename_txt = 'log_odds_matrix.txt'
    filename_csv = 'log_odds_matrix.csv'
    
    # Call the function to write to text and CSV
    log_odds_to_text_and_csv(log_score, filename_txt, filename_csv)
    
main()
