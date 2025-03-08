# Bioprolegomena

## Overview
Bioprolegomena is a university assignment that analyses nucleotide sequences using Python and additional libraries such as NumPy and Matplotlib. The script processes two nucleotide sequences and calculates key values related to sequence alignment and mutation probability.

## Features
The script performs the following calculations and visualizations:
- **Edit Distance**: Determines the best possible sequence alignment using a penalty system.
- **Mutation Probability**: Computes the likelihood of one sequence mutating into another.
- **nMax Calculation**: Identifies the highest probability under specific conditions.
- **Graphical Visualisation**: Generates a line graph showing probability values, highlighting `nMax`.

---

## Implementation Details

### 1. Interpreting Sequences from a File
This step initialises the project by extracting and validating nucleotide sequences.

#### **Process**
- Prompts the user for a text file containing two nucleotide sequences.
- Extracts sequences and assigns them to variables.
- Validates sequences for length and correctness.
- Returns valid sequences for further processing.

#### **Research & Implementation**
- Used the `linecache` library for efficient file reading.
- Incorporated `os` to check file existence and ensure it contains valid data.
- Implemented the function `sequence_checker()` to extract and validate sequences.

#### **Testing**
- Used sample sequences (`AACAGTTACC` & `TAAGGTCA`).
- Tested with incorrect inputs (invalid characters, empty files, etc.).

---

### 2. Calculating the Edit Distance
This was the most challenging part, requiring recursive and dynamic programming approaches.

#### **Process**
- Takes two sequences as input.
- Computes all possible sequence alignments using a penalty system.
- Determines the lowest penalty value (edit distance) and returns it.

#### **Research & Implementation**
- Referenced the project paper, Wikipedia’s **Sequence Alignment** article, and educational videos.
- Attempted a recursive approach but later moved to **dynamic programming** for efficiency.
- Utilised `NumPy` to optimise matrix-based calculations.

#### **Testing**
- Verified results with correct sequences.
- Tested against incorrect/malformed inputs.

---

### 3. Calculating Mutation Probability
This step determines the probability of one sequence evolving into another.

#### **Process**
- Uses mathematical models to compute mutation probability.
- Includes factors such as sequence length and mutation rate.

#### **Research & Implementation**
- Studied mutation probability models in bioinformatics literature.
- Used mathematical functions and NumPy for precision.

#### **Testing**
- Compared computed probabilities with expected values.
- Tested edge cases (identical sequences, highly different sequences).

---

### 4. Graphing Probability & nMax Calculation
The final step generates a visualisation of probability values.

#### **Process**
- Computes probabilities at different mutation steps.
- Identifies `nMax`, the highest probability point.
- Plots values using `Matplotlib`.

#### **Research & Implementation**
- Used `Matplotlib` for graphing.
- Highlighted `nMax` using annotations.
- Optimised visualisation for clarity.

#### **Testing**
- Verified graph accuracy using known probability values.
- Ensured `nMax` was correctly highlighted.

---

## Usage
1. Prepare a text file containing two nucleotide sequences.
2. Run the script and provide the file path when prompted.
3. View calculated values and generated graph.
