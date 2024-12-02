'''
This program gives us log-odd caluclations. 

This praogram takes a text file as an input of protein sequence arranged in FASTA format . It performs several analysis on amino acid sequence. 
Expected file format: 
Protein:NAME
>SequenceName1
Sequence1

>>SequenceName2
Sequence2
.
.
.

Protein:NAME2
>SequenceName1
Sequence1

>>SequenceName2
Sequence2
.
.
.

Protein:NAME3
>SequenceName1
Sequence1

>>SequenceName2
Sequence2
.
.
.

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

'''
#Used to take log 
import math
#Used to convert output to csv 
import csv
#Used to terminate program 
import sys

#Creating a list for valid amino acid including '-'
valid_amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']

"""
Function to open the FASTA format txt file to get the sequence
"""
def open_file(file_name):
    # Initializing a Dictionary to store protein families and their associated sequences
    Family_protein = {}

    #Variables to hold current family name and its sequences
    current_family = None
    sequences = []

    try:
        with open(file_name, 'r') as file:
            for line in file:
                line = line.strip()
                # If a new protein family starts, store the previous family if any
                if line.startswith('Protein:'):
                    if current_family: 
                        Family_protein[current_family] = sequences  # Saving the previous protein family
                    current_family = line.split(":")[1].strip()  # Extracting the family name after 'Protein:'
                    sequences = []  # Initializing new list for sequences
                elif line.startswith('>'):
                    continue  # Skip sequence header lines stating with > 
                elif line:  # Only add non-empty lines as sequences
                    sequences.append(line)  # Add sequence line to the current family

            # Save the last protein family if there is any remaining sequence
            if current_family:
                Family_protein[current_family] = sequences

    except IOError:
        print("Error: file does not appear to exist")
        sys.exit()

    
    return Family_protein



# Function to validate if there are any invalid characters in the sequences
def validate_aminoacid(data, valid_amino_acids):
    sequence_counter = 1  # Initializing the sequence counter to point out the number of sequence which is invalid
    for key, sequences in data.items():  # Iterating over each key-value pair in the dictionary
        for seq in sequences:  # Iterating through each sequence in the list for the current key
            if not all(char in valid_amino_acids for char in seq.upper()):  # Check the individual sequence
                print(f"Invalid amino acid sequence found in sequence number {sequence_counter} or Proteinname doesnot start with 'Protein:' ")
                exit()  # Terminating the program if an invalid sequence is found
            sequence_counter += 1  # Increment the counter for the next sequence

#QNO. 2a 
#This function calculates the frequency of single amino acid and total number of amino acid
def calculate_frequency(sequences):
    #Initializing dictionary to store the count for each individual amino acid
    counts_for_frequency = {
        'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0,
        'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0,
        'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0, '-': 0 
    }

    total_count = 0

    # Iterating through the dictionary of sequences
    for seq_list in sequences.values():
        #Second loop to go through the list inside of the dictionary
        for seq in seq_list:
            for char in seq:  # Check if the character is valid
                if char in counts_for_frequency:
                    counts_for_frequency[char] += 1  # Increment the count
                    total_count += 1  # Increment the total count
    
    # Calculate frequencies
    proportion_for_observed_frequency = {}
    # Loop through all amino acids in the dictionary
    for amino_acid, count in counts_for_frequency.items():  
        #Caculate the proportion 
        proportion_for_observed_frequency[amino_acid] = count / total_count

    return proportion_for_observed_frequency, counts_for_frequency

#Part 2b
# Function to calculate pairwise frequencies(Helper function for "calculate_amino_pair")
def calculate_pairwise(counts_for_frequency):
    # List of possible amino acids
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-']
    
    # Initialize pairwise_frequencies dictionary with all possible pairs set to 0
    pairwise_frequencies = {}
    for i in range(len(amino_acids)):
        for j in range(i, len(amino_acids)):  # Only counting (i, j) once
            pairwise_frequencies[(amino_acids[i], amino_acids[j])]  = 0 

    # Loop through the dictionary and calculate pairwise frequencies
    amino_with_count = {aa: count for aa, count in counts_for_frequency.items() if count != 0}
    for i, amino_acid1 in enumerate(amino_with_count):  # First loop (outer loop)
        count_i = amino_with_count[amino_acid1]
        
        # Now loop from `i` to the end to form pairs (i, j)
        for j, amino_acid2 in enumerate(list(amino_with_count.keys())[i:]): 
            count_j = amino_with_count[amino_acid2]
            
            # If same amino acid, calculate using the self-pair formula
            if amino_acid1 == amino_acid2:  
                pair_count = (count_i * (count_i - 1)) / 2
                # If different amino acids, calculate using the formula: count_i * count_j
            else:  
                pair_count = count_i * count_j
            
            # Accumulate the frequency into the pairwise_frequencies dictionary
            pairwise_frequencies[(amino_acid1, amino_acid2)] += pair_count

    return pairwise_frequencies

def calculate_amino_pair(sequences):
    # Initializing dictionary to store the results
    protein_results = {}
    total_pair_counts = {}  # Aggregates all unique pairs across proteins
    total_pair_count = 0  # Initializing total count for pairs

    # Iterating through each protein's sequences
    for protein_id, seq_list in sequences.items():
        # Ensuring all sequences are the same length
        length_of_individual_seq = len(seq_list[0])
        if not all(len(seq) == length_of_individual_seq for seq in seq_list):
            raise ValueError(f"Sequences for protein {protein_id} are not all the same length.")

        protein_results[protein_id] = []

        # Counting occurrences of amino acids and pairs at each position
        for i in range(length_of_individual_seq):
            # Initializing counts for this position
            counts_for_frequency = {aa: 0 for aa in 'ACDEFGHIKLMNPQRSTVWY-'}

            # Counting amino acids at position i
            for sequence in seq_list:
                amino_acid = sequence[i]
                if amino_acid in counts_for_frequency:
                    counts_for_frequency[amino_acid] += 1
            
            # Calculating pairwise frequencies for this position
            pairwise_frequencies = calculate_pairwise(counts_for_frequency)

            # Storing the result for this position
            protein_results[protein_id].append({
                "position": i,
                "pairwise_counts": pairwise_frequencies
            })

            # Aggregating the pairwise frequencies into the total counts across all positions
            for pair, count in pairwise_frequencies.items():
                if pair in total_pair_counts:
                    total_pair_counts[pair] += count
                else:
                    total_pair_counts[pair] = count
                
                # Updating the total pair count 
                total_pair_count += count  # Addinf count of the current pair to the total count
    # Add 0.5 to each pair count to avoid 0/something
    for pair, count in total_pair_counts.items():
        total_pair_counts[pair] += 0.5  

    # Loop to sum up all the values in the total_pair_counts dictionary
    total_sum = 0
    for count in total_pair_counts.values():
        total_sum += count

    # Creating the proportion list for each pair by dividing by the total_sum
    proportion_for_observed_frequency = {
        pair: count / total_sum for pair, count in total_pair_counts.items()
    }

    return protein_results, proportion_for_observed_frequency, total_pair_count

#part 2c
def expected_frequencies(aminoacid_frequency_proportion):
    # Initialize the dictionary to store expected frequencies
    expected_frequencies_dict = {}

    # Loop over all pairs of amino acids also including slefpairs
    for aminoacid1 in aminoacid_frequency_proportion:
        for aminoacid2 in aminoacid_frequency_proportion:
            # Storing the frequency (proportion) of the amino acids
            freq1 = aminoacid_frequency_proportion[aminoacid1]
            freq2 = aminoacid_frequency_proportion[aminoacid2]

            # Calculate expected frequency for the pair
            if aminoacid1 == aminoacid2:
                # For self-pairs
                expected_freq = freq1 * freq1
            else:
                # For cross-pairs 
                expected_freq = 2 * freq1 * freq2

            # Store the expected frequency in the dictionary
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

# This function outputs the final matrix in the desired format for text and CSV files
def log_odds_to_text_and_csv(log_odds_dict, filename_txt='log_odds_matrix.txt', filename_csv='log_odds_matrix.csv'):
    # Write to the text file
    with open(filename_txt, 'w') as file:
        for aa1 in sorted(set([pair[0] for pair in log_odds_dict.keys()])):  # Iterate through amino acids
            # Create a dictionary for each amino acid with the corresponding pairs and values
            row_dict = {aa2: log_odds_dict.get((aa1, aa2), 0) for aa2 in sorted(set([pair[1] for pair in log_odds_dict.keys()]))}
            # Write the formatted string in the desired structure
            file.write(f"'{aa1}': {row_dict}\n")
    
    # Write to the CSV file as usual
    # Manually create the list of unique amino acids from the dictionary keys
    amino_acids = sorted(set([aa for pair in log_odds_dict.keys() for aa in pair]))

    # Open the CSV file in write mode
    with open(filename_csv, mode='w', newline='') as file:
        writer = csv.writer(file)

        # Write the header row (amino acids as column labels)
        writer.writerow([''] + amino_acids)  

        # Write each row for each amino acid
        for aa1 in amino_acids:
            row = [aa1]  # Start the row with the amino acid label
            for aa2 in amino_acids:
                # If the pair exists in the dictionary, add the value; otherwise, add 0
                if (aa1, aa2) in log_odds_dict:
                    row.append(log_odds_dict[(aa1, aa2)])
                else:
                    row.append(0)  #0 if the pair does not exist
            writer.writerow(row)

"""
Main function to run the program logic.
"""
def main():
    #File to read
    file_name = "sequence.txt"  

    #Function to open the files and it passes sequence with name
    All_protein_sequence = open_file(file_name)

    #calling the function that validates aminoacid sequences, teminates the program if invalid amino acid sequence, or wrong protein header format
    validate_aminoacid(All_protein_sequence, valid_amino_acids)

    #Calulates the indidual and total amino acid frequency 
    #counts, individual_amino_frequencies = 
    individual_aminoacid_frequency_proportion, count_observed_frequency = calculate_frequency(All_protein_sequence)

    #calculates the amio acid pair from the sequences 
    result, proportion_for_observed_frequency, total_pair_count = calculate_amino_pair(All_protein_sequence)

    expected_frequencies_dict = expected_frequencies(individual_aminoacid_frequency_proportion)

    #final function that caclulates the log score for the final matrix 
    log_score = Logodd_scores(expected_frequencies_dict, proportion_for_observed_frequency)
 
    print(log_score)

    #Name for output files, (Change the output filename as needed)
    filename_txt = 'log_odds_matrix.txt'
    filename_csv = 'log_odds_matrix.csv'
    
    # Call the function to write into text and CSV
    log_odds_to_text_and_csv(log_score, filename_txt, filename_csv)

main()
