# COMP1730/6730 Bioinformatics Project

# YOUR ANU ID: u7637337
# YOUR NAME: Zac Morehouse

# Note that this program expects the DNA sequences to be given in a .txt document, with each sequence on line 1 and 2.
# Test cases I've been using have included in the Project.zip file. 

# Import Libraries
import os  # Used to check the file being imported
import linecache  # Used to read specific lines in a file
import numpy as np  # Used to square matrixes
import matplotlib.pyplot as mpl  # Used to plot a graph


def sequence_checker(filename, max_length):
    '''
    A function to import a text file, reading the sequences and checking they are valid.

    Args :
    filename - The name of the file being read
    max_length - The maximum length of the sequences

    Returns :
    An error message if invalid, or variables containing sequence 1 and sequence 2.
    '''

    # Base cases - check if the file given is empty or not found.
    if not os.path.exists(filename):
        raise FileNotFoundError("Error : the file does not exist.")
    elif os.stat(filename).st_size == 0:
        print("Error : the file is empty.")
        quit()

    # Ensure that our max_length is a positive integer
    if type(max_length) != int or max_length < 0:
        print("Error : max_length must be a positive integer.")
        quit()

    # Assign our sequences from line 1 and two of the given file
    sequence1 = linecache.getline(filename, 1).strip().upper()
    sequence2 = linecache.getline(filename, 2).strip().upper()

    nucleotides = ['A', 'C', 'G', 'T', '-']

    # Check the length of the sequences do not exceed the given max length
    if len(sequence1) > max_length or len(sequence2) > max_length:
        print("Error : Nucleotide sequence is too long.")
        quit()

    # Check the sequence only contains valid nucleotides
    for nucleotide in sequence1:
        if nucleotide not in nucleotides:
            print("Error : Invalid nucleotide symbol found:" + nucleotide)
            quit()
    for nucleotide in sequence2:
        if nucleotide not in nucleotides:
            print("Error : Invalid nucleotide symbol found:" + nucleotide)
            quit()

    # Check the sequences are not empty
    if len(sequence1) == 0 and len(sequence2) == 0:
        print("Error : You have not provided any sequences")
        quit()
    elif len(sequence1) == 0 and len(sequence2) != 0:
        print("Please note that Sequence 1 is empty The program will continue running, assuming S1 are all gaps.")
    elif len(sequence1) != 0 and len(sequence2) == 0:
        print("Please note that Sequence 2 is empty. The program will continue running, assuming S2 are all gaps.")

    print("Nucleotide sequences are valid and within the acceptable limit.")
    return sequence1, sequence2


def calculate_edit_distance(sequence1, sequence2):
    '''
    A function to calculate the edit distance (lowest possible penalty) when aligning two sequences, implementing the penalty system as determined by the paper.
    This function utilises a 2D table, representative of each of the nucleotides in the sequence. The function iterates over both sequences, filling the table by applying the relevant penalties (deduced by matching, inserting, deleting or substituting nucleotides to align the sequences.)
    Once the table is filled, the function takes the bottom right value (edit distance), storing it. It then passes the table to the aligner() function (with our sequences) to calculate the newly aligned sequences. 
    Finally, once the calculations are done, the function returns the aligned sequences alongside the edit distance to the user.

    Args :
    sequence1 - The first DNA sequence as determined by the sequence_checker function
    sequence2 - The second DNA sequence as determined by the sequence_checker function

    Returns :
    A tuple containing the edit distance, alongside the newly aligned sequences.
    '''

    # Create a table (list of lists), filling it initially with empty values (0). To account for empty sequences, we must + 1 to the length of the sequence
    table = [[0] * (len(sequence2) + 1) for _ in range(len(sequence1) + 1)]

    # Fill the first row and column of the table using the maximum penalty.
    for x in range(len(sequence1) + 1):
        table[x][0] = 2 * x
    for y in range(len(sequence2) + 1):
        table[0][y] = 2 * y

    # Loop through the table - going from position (1,1) to position (length of sequence 1) to (length of sequence 2).
    for x in range(1, len(sequence1) + 1):
        for y in range(1, len(sequence2) + 1):

            # First we check if the characters at each position are the same. If they are, a penalty is not incurred and the previous value is copied.
            if sequence1[x - 1] == sequence2[y - 1]:
                table[x][y] = table[x - 1][y - 1]

            # If they are not the same, we must perform three operations - substitution, deletion or insertion. We take the minimum of these values and add it to the cell as per the penalty system
            else:
                table[x][y] = min(
                    # +1 to the cell if we can substitute a nucleotide
                    table[x - 1][y - 1] + 1,
                    # +2 to the cell if we can delete a nucleotide
                    table[x - 1][y] + 2,
                    # +2 to the cell if we can insert a nucleotide
                    table[x][y - 1] + 2
                )

    # Once we have finished iterating, we pass the final table with our sequences to the aligner function - assinging them to aligned1 and aligned2
    aligned_sequence1, aligned_sequence2 = aligner(table, sequence1, sequence2)

    # Finally, return the edit distance and the aligned sequences. We do this by returning the bottom right cell in the table (minimum value)
    return table[-1][-1], aligned_sequence1, aligned_sequence2


def aligner(table, sequence1, sequence2):
    '''
    A function to calculate our newly aligned sequences by backtracking through the table filled above.
    Knowing the penalty system, we can reverse engineer the table and sequences to discover the most optimal alignment. 
    This is much less intensive than storing every possible alignment. 

    Args :
    table - The 2D table of values determined by the calculate_edit_distance function
    sequence1 - The first DNA sequence as determined by the sequence_checker function
    sequence2 - The second DNA sequence as determined by the sequence_checker function

    Returns :
    The newly aligned sequences.
    '''

    # Assign variables to the length of the sequences.
    i, j = len(sequence1), len(sequence2)

    # Create new, empty strings to hold our aligned sequences.
    aligned_sequence1, aligned_sequence2 = '', ''

    # Work backwards iterating through the sequences (iterate over until the length is 0)
    while i > 0 or j > 0:

        # If we have reached the start of sequence 1 (but not 2), add a gap, than move down to the next position
        if i == 0:
            aligned_sequence1 = '-' + aligned_sequence1
            aligned_sequence2 = sequence2[j - 1] + aligned_sequence2
            j -= 1

        # If we have reached the start of sequence 2 (but not 1), add a gap, than move down to the next position
        elif j == 0:
            aligned_sequence1 = sequence1[i - 1] + aligned_sequence1
            aligned_sequence2 = '-' + aligned_sequence2
            i -= 1

        # If both the characters match, add them to their respective sequences
        elif sequence1[i - 1] == sequence2[j - 1]:
            aligned_sequence1 = sequence1[i - 1] + aligned_sequence1
            aligned_sequence2 = sequence2[j - 1] + aligned_sequence2
            i -= 1
            j -= 1

        # Check if the current value in the table was created by substitution. If it was, we add the corresponding character to both sequences
        elif table[i][j] == table[i - 1][j - 1] + 1:
            aligned_sequence1 = sequence1[i - 1] + aligned_sequence1
            aligned_sequence2 = sequence2[j - 1] + aligned_sequence2
            i -= 1
            j -= 1

        # Check if the current value in the table was created by deletion. If it was, we add the corresponding character to seq1 and a gap to the seq2
        elif table[i][j] == table[i - 1][j] + 2:
            aligned_sequence1 = sequence1[i - 1] + aligned_sequence1
            aligned_sequence2 = '-' + aligned_sequence2
            i -= 1

        # Finally, if none of the above are true, the value must be created by insertion. If it was, we add the corresponding character to seq2 and a gap to the seq1
        else:
            aligned_sequence1 = '-' + aligned_sequence1
            aligned_sequence2 = sequence2[j - 1] + aligned_sequence2
            j -= 1

    return aligned_sequence1, aligned_sequence2


def probability(aligned_sequence1, aligned_sequence2, n):
    '''
    A function to deduce the probability of one sequence mutating to another sequence.
    This function implements a formula outlined within the paper - taking the relevant probabilities or matrix positions and multiplying them accordingly.
    This function is called many times later on to determine nmax, with n being incremented to square the matrixes. NumPy is used to simplify this process.

    Args :
    aligned_sequence1 - Aligned DNA sequence 1
    aligned_sequence2 - Aligned DNA sequence 2
    n - Number of mutations

    Returns :
    The total probability of the mutation occuring.
    '''

    # Create the matrix M using nested list.
    M = np.array([[0.976, 0.010, 0.007, 0.007],
                  [0.002, 0.983, 0.005, 0.010],
                  [0.003, 0.010, 0.979, 0.007],
                  [0.002, 0.013, 0.005, 0.979]])

    # Create dictionaries to store the relevant nucleotide probabilities and positions in the matrix.
    probability_of_bases = {'A': 0.1, 'C': 0.4, 'G': 0.2, 'T': 0.3, '-': 1}
    position_of_bases = {'A': 0, 'C': 1, 'G': 2, 'T': 3, }

    total_probability = 0

    # Multiply the matrix to the power of n as per the paper.
    new_matrix = np.linalg.matrix_power(M, n)

    # Iterate over the max length of the strings, assigning the relevant letters to variables.
    for i in range(0, max(len(aligned_sequence1), len(aligned_sequence2))):
        nucleotide1 = aligned_sequence1[i]
        nucleotide2 = aligned_sequence2[i]

        # Check if either string features a blank space (-) and substitute with one
        if nucleotide1 == '-':
            alignment_probability = 1
        elif nucleotide2 == '-':
            alignment_probability = (probability_of_bases[nucleotide1] * 1)

        # Otherwise, apply the algorithm.
        else:
            # Grab the probability of the letter, than multiply by the relevant value found in the matrix
            alignment_probability = (
                probability_of_bases[nucleotide1] * new_matrix[position_of_bases[nucleotide1]][position_of_bases[nucleotide2]])

        # If it is the first iteration, add the char_value to the total probability
        if i == 0:
            total_probability = alignment_probability
        # Otherwise, multiply it to the total probability
        else:
            total_probability = total_probability * alignment_probability

    return total_probability

def calculate_nmax(aligned_sequence1, aligned_sequence2):
    '''
    A function to calculate the nmax (the number of mutations with the highest probability)
    This function increments n, running the probability function. It then compares the results to the previous increments to deduce the highest possible probability.
    The function also stores the relevant values for graphing later on.

    Args:
    aligned_sequence1 - Aligned DNA sequence 1
    aligned_sequence2 - Aligned DNA sequence 2

    Returns :
    The highest probability of mutation, and its respective nmax value.
    '''

    nmax = 1
    n = 1
    current_probability = probability(aligned_sequence1, aligned_sequence2, n)

    # Create empty lists, appending to them for graphing later on.
    n_values = []
    likelihood_values = []

    # Keep increasing n until the current_probability stops increasing
    while True:

        # Increment n and calculate the probability.
        new_probability = probability(
            aligned_sequence1, aligned_sequence2, n+1)

        # Base cases - incase probability is always the same (empty sequences, one empty sequence, etc.)
        if current_probability == new_probability:
            print("Highest probability of mutation : 0")
            print("Maximum number of mutation steps : " + str(nmax))
            print("There is no need to graph, as n is not affecting your probability")
            quit()

        # If new probability is greater than or equal to current probability, update the current highest probability and increase n (appending the previous values)
        if new_probability >= current_probability:
            current_probability = new_probability
            n += 1

            n_values.append(n)
            likelihood_values.append(current_probability)

        # If new probability is less than the current highest probability, return the previous n value (nmax), the highest probability and break the function
        else:
        
            nmax = max(n-1, 1)

            print("Maximum number of mutation steps : " + str(nmax))
            print("Highest probability of mutation : " +
                  str(probability(aligned_sequence1, aligned_sequence2, nmax)))

            # If the probability doesn't change, don't create a graph.
            if nmax == 1:
                print(
                    "There is no need to graph, as n is not affecting your probability")
                break
            # Otherwise, create a graph.
            else:
                graph_creator(n_values, likelihood_values)
                break


def graph_creator(n_values, likelihood_values):
    '''
    A function to generate a graph using the lists defined above and MatPlotLib.

    Args:
    n_values - A list of n_values up until nmax+1
    likelihood_values - Each n_values' respective probability.

    Returns :
    A graph plotting the above.
    '''
    mpl.title('Graphing the Likelihood of Mutation against the Number of Mutation Steps')
    mpl.plot(n_values, likelihood_values)
    mpl.xlabel('n')
    mpl.ylabel('Likelihood')
    mpl.show()

# Our main function - all these lines run when the program is initialised.
if __name__ == '__main__':

    # Ask the user for the file and maximum length
    filename = input("Enter the filename containing nucleotide sequences: ")
    max_length = int(input("Enter the maximum length of the nucleotide sequences: "))

    # Read the sequences from the file
    sequence1, sequence2 = sequence_checker(filename, max_length)

    # Calculate the edit distance and aligned sequences, returning them to the user
    edit_distance, aligned_sequence1, aligned_sequence2 = calculate_edit_distance(sequence1, sequence2)
    print("Edit distance:", edit_distance)
    print("Aligned sequence 1:", aligned_sequence1)
    print("Aligned sequence 2:", aligned_sequence2)

    # Calculate the nmax with our aligned sequences. Within this, we return the nmax and probability, and graph the outcome.
    calculate_nmax(aligned_sequence1, aligned_sequence2)
