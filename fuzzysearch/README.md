## A simple console app for search for patterns or motifs in genomic sequences

Requires a path to input file with genetic sequence, a pattern and a path to output dicertory where results will be saved.
Outputs two csv files for two orientations with columns 'Sequence name', 'Start', 'End', 'Pattern searched'

Input file - FASTA or equivalent
PATTERN - The pattern string. 
Supported syntax: A, T, G, C [ATGC] or other character classes,\
                  \\D              -> any DNA base,\
                  {n}, {n,m}, {n,} -> quantifiers,\
                  |                -> alternatives.\
                  Examples:ATGC,\
                          T{1,2}A{0,1}G{3,5},\
                          T{1,2}\\D{1}A{1}G{3,5}|T{2}A{1}G{2}')


```
usage: script.py -f <path> -s <pattern> -d <dir>

Search for motif or pattern and save results.

options:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  Path to input file
  -s PATTERN, --pattern PATTERN
                        The pattern string

  -d DIR, --dir DIR     Path to results directory

```
