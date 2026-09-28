# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "biopython",
#     "pandas",
# ]
# ///
import argparse
import re
import pandas as pd
import csv
import os
from Bio.SeqIO.FastaIO import FastaIterator
from pathlib import Path

"""
usage: script.py -f <path> -s <pattern> -d <dir>

Search for genomic or text patterns and save results.

options:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  Path to input file
  -s PATTERN, --pattern PATTERN
                        The pattern string

  -d DIR, --dir DIR     Path to results directory

"""


parser = argparse.ArgumentParser(description="Process a file and a pattern.")

# Add arguments
parser.add_argument('-f', '--file', type=str, required=True, help='Path to file')
parser.add_argument('-s', '--pattern', type=str, required=True, help=
                    'The pattern string. Supported syntax:\
                                A, T, G, C\
                                [ATGC] or other character classes,\
                                \\D              -> any DNA base,\
                                {n}, {n,m}, {n,} -> quantifiers,\
                                |                -> alternatives.\
                            Examples:\
                                ATGC,\
                                T{1,2}A{0,1}G{3,5},\
                                T{1,2}\\D{1}A{1}G{3,5}|T{2}A{1}G{2}')

parser.add_argument('-d', '--dir', type=str, required=True, help='Path to results')

args = parser.parse_args()

# Access the variables

path_object = Path(args.file)
print(f"Input file exists? {path_object.exists()}")

pattern_object = args.pattern
if not pattern_object.strip():
    print("Error: The pattern cannot be empty or just spaces.")
else:
    print(f"Pattern successfully received: {pattern_object}")

output_object = Path(args.dir)
print(f"Output directory? {output_object.exists()}")



def reverse_complement_regex(pattern):
    """
    Reverse-complement a DNA motif regex.
    """

    complement_map = str.maketrans("ATGC", "TACG")

    token_pattern = re.compile(
        r"""
        \[[^\]]+\](?:\{[^}]+\})?   # character class
        |
        \\D(?:\{[^}]+\})?           # DNA wildcard
        |
        [ATGC](?:\{[^}]+\})?        # nucleotide
        """,
        re.VERBOSE
    )

    def reverse_branch(branch):
        tokens = token_pattern.findall(branch)
        return "".join(
            token.translate(complement_map)
            for token in reversed(tokens)
        )

    return "|".join(
        reverse_branch(branch)
        for branch in pattern.split("|")
    )

patternG = str(pattern_object).strip().replace(' ', '').upper().replace('U', 'T')
patternC = reverse_complement_regex(patternG)

patternG_obj = re.compile(patternG)
patternC_obj = re.compile(patternC)

#Parsing fasta 
def fasta_reader(infile):
    with open(infile) as handle:
        for record in FastaIterator(handle):
            yield record

#Search pattern
def search_functionGC(myseq, pattern):
    ind_spans = [match.span() for match in re.finditer(pattern, myseq)]

    ind_start=[span[0] for span in ind_spans]
    ind_end=[span[1] for span in ind_spans]

    return ind_start, ind_end


def repeats_dataframeGC(id_, ind_start, ind_end):

    pattern_df = pd.DataFrame({   
        "Sequence name": [id_]*len(ind_start),
        "Sequence Start": ind_start,
        "Sequence End" : ind_end})

    return pattern_df

def csvGC_repeats_open(patternG, patternC, save_dir, path_object):
    file_stem = path_object.stem 
    
    safe_patternG_name = patternG.replace('|', '_')
    safe_patternC_name = patternC.replace('|', '_')
    invalid_chars = ['/', '\\', ':', '*', '?', '"', '<', '>', '[', ']']
    for char in invalid_chars:
        safe_patternG_name = safe_patternG_name.replace(char, '')
    for char in invalid_chars:
        safe_patternC_name = safe_patternC_name.replace(char, '')

    csv_path_repeatsG = os.path.join(save_dir, f'repeats_{file_stem}_{safe_patternG_name}.csv')
    csv_file_repeatsG = open(csv_path_repeatsG, 'w', newline='')
    headersG_repeats = ['Sequence name', 'Start', 'End', f'Pattern searched {patternG}']
    writer = csv.writer(csv_file_repeatsG, delimiter=',')
    writer.writerow(headersG_repeats)
    
    csv_path_repeatsC = os.path.join(save_dir, f'repeats_{file_stem}_{safe_patternC_name}.csv')
    csv_file_repeatsC = open(csv_path_repeatsC, 'w', newline='')
    headersC_repeats = ['Sequence name', 'Start', 'End', f'Pattern searched {patternC}']
    writer = csv.writer(csv_file_repeatsC, delimiter=',')
    writer.writerow(headersC_repeats)

    return csv_path_repeatsG, csv_path_repeatsC

def save_repeats_csvGC(repeatsG, repeatsC, csv_file_repeatsG, csv_file_repeatsC):
    repeatsG.to_csv(csv_file_repeatsG, mode='a', sep=',', index=False, header=False, date_format=object)
    repeatsC.to_csv(csv_file_repeatsC, mode='a', sep=',', index=False, header=False, date_format=object)

    return repeatsG, repeatsC

repeats_files = csvGC_repeats_open(patternG, patternC, output_object, path_object)
csv_file_repeatsG = repeats_files[0]
csv_file_repeatsC = repeats_files[1]

for entry in fasta_reader(path_object):
    seqname = entry.id
    seq = entry.seq
    description = entry.description
    print(description)
    myseq = str(seq).replace(' ', '').upper().replace('U', 'T')
    resultG = search_functionGC(myseq, patternG_obj)
    resultC = search_functionGC(myseq, patternC_obj)

    ind_startG = resultG[0]
    ind_endG = resultG[1]
    ind_startC = resultC[0]
    ind_endC = resultC[1]


    repeatsG = repeats_dataframeGC(seqname, ind_startG, ind_endG)
    repeatsC = repeats_dataframeGC(seqname, ind_startC, ind_endC)

    saved_repeats = save_repeats_csvGC(repeatsG, repeatsC, csv_file_repeatsG, csv_file_repeatsC)
