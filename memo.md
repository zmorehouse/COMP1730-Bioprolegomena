# Bioinformatics Project Memo
**Zac Morehouse**
**Student Number : u7637337**

## Brief Intro
As someone coming from a background in Humanities - this project was originally very daunting. Not only was I unfamiliar with general Biology and Physics concepts (Projects 3 and 2), but I also struggled to understand the higher Mathematical concepts outlined in Project 1. Consequently, I opted for Paper 3, with the project outline suggesting that it would be the best choice for someone with no prior knowledge in the fields. 

Though Project 3 was still quite challenging, I found it to be much more manageable when broken down into smaller problems and collated together at the end. Resultingly, this report is structured as such - outlining each 'peice' of the problem, how they were researched, implemented and tested before being collated together to produce my final solution.

## Breaking the Project Down
After initially reading through the project, I split the project into the following peices :

1. Make the program read and interpret two given Nucleotide sequences.
2. Calculate the edit distance (best possible sequence alignment) and return to the user.
3. Calculate the two best aligned sequences and return to the user.
4. Calculate the probability of one string mutating to another.
5. Calculate the nmax(highest probability) that can be obtained and return to the user.
6. Keep track of steps 4 and 5. Graph these using matplotlib.

## 1. Interpreting Sequences from a File.
This peice of the project required some initial setup for use further on in the project. It was more-or-less simple enough, though researching external libraries made this process much smoother.

**What this step needed to do**
- Prompt the user for a text file containing 2 nucleotide sequences, alongside the maximum length of said sequences.
- Prompt the user for the maximum length of the strings
- Extract each sequence individually and assign to variables
- Check them to ensure they were the correct length, had the correct values, etc.
- If correct, return the sequences as variables.

**Researching**
Little external research was neede here - the paper clearly outlined the requirements, and additional base cases could be inferred from this. After a quick look online, I decided implementing the [linecache library](https://docs.python.org/3/library/linecache.html) would be useful - both reducing difficulty of implementing the code, and making it easier to read. I also added the [os library](https://docs.python.org/3/library/os.html) to check if the file could be found, and that it had values in it.

**Implementing** 
Following the above, I defined the function `sequence_checker`. Within this, I used the os and linecache libraries to read the given file, checked the sequences against some intial criteria (its contents, length, etc.) and returned the sequences as variables for the next step in the program.

**Testing** 
- Tested with a variety of given examples (AACAGTTACC & TAAGGTCA)
- Tested with intentionally broken examples (incorrect letters, empty files, etc.) 

## 2. Calculating the Edit Distance 
By far the most difficult part of the assignment. Not only did I struggle alot following the recursive method, but attempting to implement a dynamic programming approach was also quite difficult.

**What this step needed to do**
- Take two sequences as arguments and compare them
- Using the penalty system, calculate the penalties of creating all possible sequence alignments. 
- From these calculations, take the lowest penalty value (edit distance) and return it to the user.

**Researching**
Initially, I read through the paper given with the project - additionally looking at the [Wikipedia page for Sequence Alignment](https://en.wikipedia.org/wiki/Sequence_alignment), alongside a video on [Sequence Alignment](https://www.youtube.com/watch?v=9bCkAsaP_z4) to understand the context of the problem. Once I had an understanding of what was needed, I begun working on a recursive approach. To help, looked to the [W3Schools Python Function Recursion Page](https://www.w3schools.com/python/gloss_python_function_recursion.asp), alongside an educational Youtube Video on [Recursive Programming](https://www.youtube.com/watch?v=ixdr6V2vRC4). 

After being unsuccessful with my first attempt at a solution, I looked further into [Memoization](https://stackoverflow.com/questions/73536388/dynamic-programming-smallest-cost-path-through-matrix-memoization), before stumbling onto a [Paper from MIT discussing Dynamic Programming in relation to Sequence Alignment (notably, pages 12 - 16)](http://web.mit.edu/6.047/book-2012/Lecture02_DynamicProgramming/Lecture02_DynamicProgramming_standalone.pdf). Though the paper didnt feature any code outright (and implemented a different penalty system), it did outline how a matrix could be used to calculate and store penalties (in relation to DNA sequencing.) This was extremely useful, as it became the basis to the logic that my dynamic approach was created on. This method proved sucessful, being implemented in my final solution.

**Implementing**
At first, I tried to implement a recursive function (`calculate_edit_distance`) that would iterate over the sequences, running each alignment technique (substitution, deletion and insertion) and taking the minimum value. A snippet of what I was attempting to acheive can be seen below. 

```python 
    # Initial attempts at a recursive function
    substitution, s1_sub, s2_sub = calculate_edit_distance(sequence1[1:], sequence2[1:])
    substitution += 0 if sequence1[0] == sequence2[0] else 1
 
    insertion, s1_ins, s2_ins = calculate_edit_distance(sequence1, sequence2[1:])
    insertion += 2

    deletion, s1_del, s2_del = calculate_edit_distance(sequence1[1:], sequence2)
    deletion += 2

    if substitution <= insertion and substitution <= deletion:
        return substitution, sequence1[0] + s1_sub, sequence2[0] + s2_sub
    elif insertion <= deletion:
        return insertion, '-' + s1_ins, sequence2[0] + s2_ins
    else:
        return deletion, sequence1[0] + s1_del, '-' + s2_del
```
Evidently, approaching the problem with a recursive methodology proved to be quite difficult. I was often unable to keep track of where my program was at and had many difficulties troubleshooting and bugfixing.

After giving up on the above, and researching online, I decided to give dynamic programming a crack - initialising 2D nested lists and storing the outcomes of the various alignment penalties within it. `table = [[0] * (len(sequence2) + 1) for _ in range(len(sequence1) + 1)]` was used to generate an empty 'table' of nested lists, which is then iterated over - filling the 'cells' with each respective penalty. The usefulness of having the numbers being stored meant that we could reference a previous cell (often the value at `[i - 1][j - 1]`) in our calculations, drastically reducing computational load. 

Once the table is filled, we can find our edit distance in the bottom right of the table `[-1][-1]` Overally, this approach - though still difficult - was much more manageable than the previous.

**Testing**
After implementing the code, I tested using a handful of examples as outlined in the paper. These included : AACAGTTACC and TAAGGTCA, making sure that 7 was being returned - but that 8 was also featured in an alternate alignment. I also tested some basecases, including -GTAA and AG--C to ensure gaps were being handled correctly and that gaps at the beginning of sequences would not break the result.

## 3. Calculating the Aligned Sequences
Not **as** tricky as the above, but still tested my understanding of the implemented algorithms and data structures.

**What this step needed to do**
- Take the data previously calculated
- Backtrack through the previous data, using its most ideal path to work out the aligned sequences
- Work out when the sequences were full and return them to the user.

**Researching**
Thankfully, all the external research from the above applied to this peice of the program (with the [Paper from MIT](http://web.mit.edu/6.047/book-2012/Lecture02_DynamicProgramming/Lecture02_DynamicProgramming_standalone.pdf) again coming in useful.)

**Implementing**
Originally I had attempted to store the various sequences in nested lists as the `calculate_edit_distance` function ran - although this proved to be far too computationally intensive. Moreover, there was the issue of keeping track of each sequence and aligning it to the corresponding penalty. 

In the end (after implementing the dynamic approach), I decided to back-track through the matrix of values created earlier. The nature of the structure means that the edit distance at any given time runs on the diagonal of the table. Using this logic, I was able to implement a new function, `aligner` which would iterate over the length of the initial sequences, finding each diagonal coordinate (`table[x][y]`), comparing it to the relevant sequence and adding it's value to its respective strings `aligned_sequence1` and `aligned_sequence2`. 

Again, akin to the above, the following cases were possible : Either the nucleotides **matched**, the value was created by **substitution**, the value was created by **deletion** or value was created by **insertion**. As such, the function features a large if statement to address these cases.

Notably, when I first wrote the function - I ran into a number of issues, as I had failed to account for differently lengthed strings. Resultingly, I also added the below code to check if we had reached the start of the strings, adding the relevant gaps.

```python
if i == 0:
    aligned_sequence1 = '-' + aligned_sequence1
    aligned_sequence2 = sequence2[j - 1] + aligned_sequence2
    j -= 1 

elif j == 0:
    aligned_sequence1 = sequence1[i - 1] + aligned_sequence1
    aligned_sequence2 = '-' + aligned_sequence2
    i -= 1
```

**Testing**
Here, I looked to test the program with the all the cases mentioned in the step above. These included : AACAGTTACC & TAAGGTCA (as outlined in the paper), CGTAA & AGTAC, ATGCATTAAGCCAGTCA & ATGCGTAAGCACTCA and others. Some base cases were also run to test gaps in sequences and differently lengthed sequences.

## 4. Calculating the probability of one string mutating to another.
I didn't struggle here as much as I did the edit distance. In fact, I actually opted to tackle this part of the assignment **FIRST** using the test cases outlined in the document. 

**What this step needed to do**
- Store the values outlined in the paper (Matrix, Nucleotide Probability, Nucleotide Position)
- Iterate over the nucleotides in each sequence
- Apply the formula M~S1 Nucleotide Probability~ * M~[S1 Nucleotide Position][S2 Nucleotide Position]~
- Add the result to a tracked total.

**Researching** 
I first looked into [NumPy](https://numpy.org/) as outlined in the paper - however, didn't see much need for it. This would come back to bite me however, as I later discovered the matrixes would need to be placed to the power of (more on that below). Outside of this, little external research was required.

**Implementing**
Firstly, I created the Matrix as defined within the paper. This consisted of nesting 4x lists within a larger list, and placing the respective values at the desired positions. 
From there, each nucleotide needed to be assigned two relevant values - a probability value and a positional value (related to the matrix). Initially, I tried to implement these using tuples however, as my code progressed, I found using dictionaries / key value pairs alot friendlier. Resultingly, I created `probability_of_bases` and `position_of_bases` dictionaries to reference from in my calculations.
Once these data sets had been created it was as simple as defining the function `probability`, looping over the length of the sequences and calculating the relevant formulas as outlined by the paper.

**Testing**
- Tested using the two given sequences in the paper : CCAT & CCGT = 0.000300 | CCGAT & CC–GT = 0.0000600

## 5. Calculating the nmax
This was a little more difficult - however, having the above `probability` function already written made it alot more manageable.

**What this step needed to do**
- Calculate the Matrix to the power of n
- Run our aligned sequences with the predetermined algorithm against a variety of incrementing Matrixes (rather than just one)
- Check the results to see which n gives us the largest probability of mutation (nmax)
- Return both the highest probability and the nmax to the user.

**Researching** 
As mentioned above, after reading further into it, I discovered I'd need to be multiplying matrixes. Initially, I thought this was as simple as multiplying the numbers by themselves, however, after [further research on Khan Academy](https://www.khanacademy.org/math/precalculus/x9e81a4f98389efdf:matrices/x9e81a4f98389efdf:multiplying-matrices-by-matrices/v/matrix-multiplication-intro) - I discovered that matrix multiplication was entirely different to your standard multiplication (and alot more difficult). 
Luckily, after looking again into [NumPy](https://numpy.org/doc/stable/reference/routines.linalg.html), I discovered that the library had a built in function `linalg.matrix_power()` that would return the square of a matrix. Resultingly, the NumPy was imported, with the pre-existing matrix being adjusted. 

**Implementing**
Firstly, I imported the NumPy Library adding `np.` to the matrix set up above as to call it in the future.
From there, I started by modifying the `probability` function, passing a new argument (`n`)  to determine what to square the matrix to. Using `new_matrix = np.linalg.matrix_power(M, n)` I was able to assign the original matrix (M) to new_matrix, placing it to the power of n to use in our calculations.
After this, I created a new function `calculate_nmax` which took our two sequences and ran them through the `probability` function. This function looped - taking the result of `probability`, and comparing it to the previously calculated result, incrementing n by one each time.
Using an if statement, the function then checks whether the newly calculated result is less than the previous. If it is, this means the previous result (`n-1`) is our nmax. 
Finally, I added a print statement to return both the calculated nmax, alongside our highest probability to the user.

**Testing**
Testing here proved to be quite difficult. As the paper didn't give any examples of nmaxes, I had to calculate these myself. Resultingly, a number of times I was unsure whether it was my program - or rather my maths - that was off. Eventually however, using a small handful of the test cases above, my program seemed to product the correct results.

## 6. Graphing the result
Thankfully, I was already somewhat familiar with matplotlib from previous lab exercises. Resultingly, this was more-or-less straightforward.

**What this step needed to do**
- Keep track of the generated n's and their respective probability values
- Plot these values against one another, with n being the x axis and probability being the y axis.
- Return a graph to the user with this information.

**Researching**
The main thing that needed to be researched was our required mpl functions. A quick look on the [MatPlotLib Documentation](https://matplotlib.org/stable/index.html) and we were good to go. No other external research was required here.

**Implementing**
Firstly, the program needed to store all the data points, rather than only returning the final result. To do this, I created two empty lists (`n_values = []` and `likelihood_values = []`).
From there, I appended them with their respective values each time the `probability` function was run.
After this, I defined a new function `graph_creator`, which took each of the lists. Calling Using `plot(n_values, likelihood_values)`, I was able to plot the values stored in the list. I also added some additional labels for context sake, before calling the function from within the `calculate_nmax` function.

**Testing**
- Added some dummy variables to the lists and ensured they were appearing correctly
- Using the VSCode debugger, checked the list variables were accurate
- Tested with : ATGCATTAAGCCAGTCA & ATGCGTAAGCACTCA alongside AACAGTTACC & TAAGGTCA

## Final Steps
Once I had peiced all of the above together, I ran some tests to ensure all was working as intended. After that, I looked to clean some of my code up by doing the following :

1. Looked through and removed any redundant code that may be present
2. Branched some additional functions to make the code easier to understand
3. Renamed some variables and functions to ensure consistency through the program
4. Added a __main__ conditional block to again make the code easier to understand
5. Added some additional comments and relevant docstrings to all functions.

And that was the end of it. Overally, this project was quite challenging, but definently helped improve my understanding of data structures, recursive and dynamic programming and the methodology through which programs should be simplified and tackled.