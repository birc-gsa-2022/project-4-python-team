[![Open in Visual Studio Code](https://classroom.github.com/assets/open-in-vscode-c66648af7eb3fe8bc4f294546bfd86ef473780cde1dea487d3c4ff354943c9ae.svg)](https://classroom.github.com/online_ide?assignment_repo_id=9324980&assignment_repo_type=AssignmentRepo)
# Project 4: BWT-based matching (FM-index)

Now that you have a functioning suffix array, you should implement the BWT-based search, also known as FM-index. This algorithm improves the search in the suffix array from O(m log n + z) to O(m + z) after O(n) preprocessing (plus whatever time it takes you to build your suffix array).

You should implement a suffix array construction algorithm. You can choose to implement the naive algorithm where you explicitly sort strings, or the O(n) skew or SAIS algorithms, or any other algorithm. After constructing the suffix array, you should implement the binary search based and the Burrows-Wheeler based search algorithm.

The algorithms should be implemented in a program named `fm`. Since we are building data structures in a preprocessing step, and since a common usage of read mappers is to map multiple files of reads against the same genome, we should build the tool such that we can preprocess a genome once, and then reuse the preprocessed data on subsequent searches.

Therefore, your tool should have options for either preprocessing or read-mapping. If you run it as `fm -p genome.fa` it should preprocess the sequences in `genome.fa`, and if you run the tool as  `fm genome.fa reads.fq` it should search the genome and produce output in the same format as the previous projects.

When you preprocess `genome.fa` you should write the result to file. You are free to choose what you write to file, how many files you use, or how you represent the output. Use the input file name, here `genome.fa` but it can be any file name, to select the file names for your preprocessed data. That way, when you run a search with `fm genome.fa reads.fq`, your tool can determine which preprocessed files to read from the second first argument.

## Evaluation

Once you have implemented the `fm` program (and tested it to the best of your abilities) fill out the report below, and notify me that your pull request is ready for review.

# Report

## Preprocessing

*What preprocessing data do you store in files, and how?*
For pre-processing we ended up using pickle, there was some amount of fidling around with storing stuff in simple txt files. But hey, why invent a way to pack and parse stuff when those already exist. 
We store all the necessary information in a class called fm index, this includes the Suffix array, the string (as a rotating string), a mapping between alphabet & index, and the O and C tables. 
The pre processing calculates these and dumps them to a binary file using pickle.

The data necessary data structures should follow these complexities:
* Suffix array - O(n log n) no, we don't have a linear implementation yet...
* Rotating string - O(n), unfortunately we make a copy of the input string, but hey this was easier
* O table - O(n) given that we can assume a constant sized alphabet, as the table is sized n * $\sigma$ 
* C table - O(n) we do one pass over the alphabet, one pass over the string for counts. we then celebrate passover.
* Alphabet mapping - O(n) one iteration on the alphabet...

That linear time SA construction is looking very appealing right about now...




## Insights you may have had while implementing the algorithm
We played around abit with RLE and the other tricks related to the bwt, so if we wanted to reduce the size of our pre-processed data, we could drop the string, and replace it with the bwt l-column compressed with rle, which would save us a bit. but that would then cost in loading time. On the same note we could cut out the suffix array, as that can be reconstructed in reasonable time (we atleast convinced ourselves that it should be possible in O(n) using the O table). It looks like there is ample opportunity for playing around with compression.

## Problems encountered if any
We didn't encounter any soulcrushing problems, only minor soulcrushing hiccups. 
Jokes aside, this proved to be one of the easier implementations we've done, so as far as we're aware we didnt encounter any major problems

## Validation

*How did you validate that the preprocessing and the search algorithm works?*
Once again, for searching we can compare results with our earlier implementations. 
Generally while testing if the algorithm worked we ensured that it would handle empty strings as both x and pattern, aswell as single letter, and repetetive strings. So far none of these break.

The pre-processing was tested by dumping an object to file and reloading, then comparing with a newly constructed to see if anything went missng, and by running pre-processing + search in the same vs in seperate scripts.

## Running time

*List experiments and results that show that both the preprocessing algorithm and the search algorithm works in the expected running time. Add figures by embedding them here, as you learned how to do in project 1.*
