# GenomicSequenceAnalyzer

This is a prototype version. Future updates will provide siginificantly more customization options for enhanced data analysis.

---------------------------------------------------------------------------------------------------------------------------------------------------

# Program Overview - Purpose and Functions

This program is an analysis tool that performs several functions to aid researchers who are working with large sequence datasets.
It divides sequence data into labeled segments and categorizes them by feature type according to the annotations in the input file.
It is intended to highlight specific areas of interest within entire genomes by filtering out insignificant regions.

Its main feature is sorting segments by a desired measure(ie. order of appearance) and then performing range queries on this sorted data.
Its flexibility with sorting measures makes this useful for any researcher as they can adjust the sorting to accommodate for their research.

---------------------------------------------------------------------------------------------------------------------------------------------------

# Data Structure Implementations 

This tool uses binary search trees and hashmaps as its primary data structures. Binary search trees are used to store sorted (and possibly 
filtered) data. We use hashmaps to store binary search trees with different filtered segments and to store the contents of the original sequence.

This program utilizes the Biopython package to parse the GenBank file's sequence and annotations. The main annotations we retrieve are the label, 
feature type, and the start and end positions. This parsing is coupled with an algorithm that stores each segment in a tree node with its details
and inserts the node into a BST by a specified measure.

This program implements a self-balancing binary search tree, or AVL tree. This ensures that no matter what order the segments are inserted into 
the tree or what sorting measure the user selects, the end result will be a relatively balanced tree, allowing for optimized range query efficiency.

---------------------------------------------------------------------------------------------------------------------------------------------------

# How to Use

Required Files:

sequenceAVLTree.c - Contains logic for AVL tree and range query operations.

filterSequence.c - Contains logic for filtering segments.

sequenceAVLTree.h - Header file containing data structure definitions and function declarations.

parse_genbank.py - A python script that parses GenBank files and creates a csv file with details for the C program to use.

Makefile - Used to compile the C source files.

.gb file - The input genomic sequence data.
  
Instructions:

  - download all required files listed above
      - select a .gb file from the provided input folder or
      - download a .gb file of interest from GenBank

  - install Biopthon (for parse_genbank.py)
      - try 'pip install biopython' or 'pip3 install biopython' on your command line
   
  - run 'make' on the command line
  - to start the program, run './sequenceAVLTree mySequence.gb' on the commmand line, replacing 'mySequence' with the name of your .gb file

  - Select a command from the list to start analyzing your sequence

---------------------------------------------------------------------------------------------------------------------------------------------------

# Current Limitations and Future Implementation Plans

Input file

Current Limitations:     only accepts .gb files(GenBank) for sequence data because of its annotation format(labels and feature types)

Future Implementations:  accept various file types that contain similar or more details

Sorting

Current Limitations:     only sorts labeled segments by order of appearance

Future Implementations:  planning to implement sorting by segment length, nucleotide/amino acid composition, pattern frequency..

Duplicate Segments

Current Limitations:     duplicate segments with different feature types excluded from trees but included in the main sequence hash

Future Implementations:  avl trees will store all feature types of a unique segment, not just the first appearance

Filtering

Current Limitations:     only filters out segments by start positions

Future Implementations:  filter out by segment labels, feature types, positional ranges..

Range Query

Current Limitations:     only accepts start positions as upper and lower bounds

Future Implementations:  accept segment labels or any sequence positions as upper and lower bounds

User Interface

Current Limitations:     not much flexibility with input types for commands

Future Implementations:  option to export filtered data or range query results to output text file

---------------------------------------------------------------------------------------------------------------------------------------------------
