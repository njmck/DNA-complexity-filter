# DNA Complexity Filter

A group of Python functions that can remove complexities of DNA sequences for submission to IDT gBlocks or Genscript Gene Synthesis, who have strict requirements for GC content and repetitive DNA sequences. Detection of DNA hairpins and palindromic sequences were planned but are not a priority currently.

In even just a 201bp DNA sequence, the number of permutations in terms of codon combinations is astronomically high. This script will randomly generate different codon combinations and accept or reject them based on the arguments you provide.

## DNA_complexity_filter.py functions:

### DNA_to_protein(DNA_seq)
1. Accepts a DNA sequence string.
2. Returns the translated protein sequence as a list of three-letter amino acids.

### GC_perc(DNA_seq)
1. Calculate the GC content of a DNA sequence string.
2. Accepts a DNA sequence string.
3. Returns the GC percentage of the DNA sequence as a float.

### DNA_filt(DNA_seq):
1. Filters the DNA string for non-ATGC characters and for divisibility by 3.
2. Returns a filtered DNA string with U changed to T.

### random_codon_assembly(aa_seq):
1. Accepts a sequence of amino acids.
2. Generates a DNA sequence with randomly-selected codons corresponding to the amino acid.

### DNA_repeats_filter(seq, max_repeat_length)
1. Accepts a DNA sequence string and scans for repeats based on the maximum repeat length provided.
2. Returns True or False depending if the sequence passes the test parameters.

### GC_sections(DNA_seq, section_len, GC_section_lower, GC_section_upper)
1. Filters a DNA sequence for GC content in each stretch of n nucleotides.
2. Returns True or False depending whether the test passes.

### filtered_codons(seq, no_of_results, max_rand_seq, GC_lower, GC_upper, GC_section_len, GC_section_lower, GC_section_upper, max_repeat_length):
* The most useful function that depends on many of the above functions.
* As a starting point, use the following values, and adjust as necessary: filtered_codons(seq, 10, 1000000, 0, 70, 100, 0, 75, 24)
* seq = A DNA sequence or list containing three-letter abbreviations of amino acids
* no_of_results = number of desired results to be returned in the list of filtered DNA sequences. If max_rand_sequences is reached first, you won't get your desired number of sequences.
* max_rand_seq = number of randomly-generated DNA sequences to generate. For tough sequences, the higher the better. But, this is time consuming and CPU intensive. If no_of_results is satisfied first, max_rand_seq will not be reached and the script stops.
* GC_lower = The lower threshold for GC% across the entire DNA sequence.
* GC_upper = The upper threshold for GC% across the entire DNA sequence.
* GC_section_len = The length of the sub-sequences within the query DNA sequence that will be tested separately for GC content.
* GC_section_lower = The lower threshold for GC% across n number of bases within the sub-sequence of specified length.
* GC_section_upper = The upper threshold for GC% across n number of bases within the sub-sequence of specified length.
* max_repeat_length = The maximum length of a repetitive sequences that are allowed to appear across the entire DNA sequence.
1. Accepts a list of three-letter residues or DNA sequence.
2. Returns a list of resulting sequences with filtered GC-content and repetitive sequences.