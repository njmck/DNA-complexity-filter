import copy
import random
import sys

# MAIN PIPELINE STARTS HERE:

DNA_string = "CGGGGATTCTCTTGAATTTATTGCATGGCCGCTGCAATTACCGATATGGCTGATCTCGAGGAGTTGTCCAGGCTCAGCCCATTGCCCCCCGGTAGTCCGGGCAGTGCTGCCAGAGGTAGAGCTGAGCCACCAGCAGCTGCTGCAGCTGCTGCAGCAGCAGCCGCAGCTGCAGAGGCAGAAGCAGTGGCAGCCCTGCTTTTG"

codon_dict = {"Phe": ["TTT", "TTC"],
              "Leu": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
              "Ile": ["ATT", "ATC", "ATA"],
              "Met": ["ATG"],
              "Val": ["GTT", "GTC", "GTA", "GTG"],
              "Ser": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
              "Pro": ["CCT", "CCC", "CCA", "CCG"],
              "Thr": ["ACT", "ACC", "ACA", "ACG"],
              "Ala": ["GCT", "GCC", "GCA", "GCG"],
              "Tyr": ["TAT", "TAC"],
              "His": ["CAT", "CAC"],
              "Gln": ["CAA", "CAG"],
              "Asn": ["AAT", "AAC"],
              "Lys": ["AAA", "AAG"],
              "Asp": ["GAT", "GAC"],
              "Glu": ["GAA", "GAG"],
              "Trp": ["TGG"],
              "Cys": ["TCT", "TGC"],
              "Arg": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
              "Gly": ["GGT", "GGC", "GGA", "GGG"]}

def DNA_to_protein(DNA_seq):
    '''
    Accepts a DNA sequence string.
    Returns the translated protein sequence as a list of three-letter amino acids.
    '''
    prot_sequence = []
    codon_sequence = []
    aa_length = len(DNA_seq) / 3
    GC_str_copy = copy.deepcopy(DNA_seq)
    aa_length_countdown = copy.deepcopy(aa_length)
    while aa_length_countdown > 0:
        GC_string_codon = GC_str_copy[:3]
        for i_0 in codon_dict:
            for i_1 in codon_dict[i_0]:
                if GC_string_codon == i_1:
                    prot_sequence.append(i_0)
                    codon_sequence.append(i_1)
                    GC_str_copy = GC_str_copy[3:]
                    aa_length_countdown -= 1
    return(prot_sequence)

def GC_perc(DNA_seq):
    '''
    Calculate the GC content of a DNA sequence string.
    Accepts a DNA sequence string.
    Returns the GC percentage of the DNA sequence as a float.
    '''
    if isinstance(DNA_seq, str):
        G_count = DNA_seq.count("G")
        C_count = DNA_seq.count("C")
        GC_count = G_count + C_count
        percent = (GC_count / len(DNA_seq)) * 100
        return (float(percent))
    else:
        print("Error. Provided argument is not type 'string'.")

def DNA_filt(DNA_seq):
    '''
    Filters the DNA string for non-ATGC characters and for divisibility by 3.
    Returns a filtered DNA string with U changed to T.
    '''
    if ((len(DNA_seq) % 3) != 0):
        print("Incompatible Sequence. Indivisible by 3. Edit your DNA sequence and try again.")
        return(None)
    filtered_DNA_seq = DNA_seq
    if "U" in list(DNA_seq):
        DNA_list = list(DNA_seq)
        for n_0, i_0 in enumerate(DNA_list):
            if i_0 == "U":
                DNA_list[n_0] = "T"
        filtered_DNA_seq = ("").join(DNA_list)
    DNA_list = list(filtered_DNA_seq)
    accepted_nucl = ["A", "T", "G", "C"]
    for i_0 in DNA_list:
        if i_0 not in accepted_nucl:
            print("Incompatible Sequence. Sequence must consist of A, T, U, G, and C.")
            return(None)
    return(filtered_DNA_seq)

def random_codon_assembly(aa_seq):
    '''
    Accept a sequence of amino acids.
    Generates a DNA sequence with randomly-selected codons.
    '''
    rand_codon_list = []
    for i_0 in aa_seq:
        rand_codon = random.choice(codon_dict[i_0])
        rand_codon_list.append(rand_codon)
    rand_codon_list_join = ("").join(rand_codon_list)
    return(rand_codon_list_join)

def DNA_repeats_filter(seq, max_repeat_length):
    '''
    Accepts a DNA sequence string and scans for repeats based on the
    maximum repeat length provided.
    Returns True or False depending if the sequence passes the test parameters.
    '''
    snippit_list = []
    bottom_num = 0
    top_num = max_repeat_length - 1
    repeat_test_pass = True
    while top_num != (len(seq) - 1):
        snippit = seq[bottom_num:top_num]
        if snippit not in snippit_list:
            snippit_list.append(snippit)
            bottom_num += 1
            top_num += 1
        else:
            repeat_test_pass = False
            break
    return(repeat_test_pass)

def GC_sections(DNA_seq, section_len, GC_section_lower, GC_section_upper):
    '''
    Filters a DNA sequence for GC content in each stretch of n nucleotides.
    Returns True or False depending whether the test passes.
    '''
    arg_length = len(DNA_seq)
    if section_len > arg_length:
        if GC_section_lower < GC_perc(DNA_seq) < GC_section_upper:
            return(True)
    else:
        bottom_num = 0
        top_num = copy.deepcopy(section_len)
        while top_num != (len(DNA_seq) - 1):
            snippit = DNA_seq[bottom_num:top_num]
            if GC_section_lower < GC_perc(snippit) < GC_section_upper:
                bottom_num += 1
                top_num += 1
            else:
                return(False)
        return(True)

def filtered_codons(seq, no_of_results, max_rand_seq, GC_lower, GC_upper, GC_section_len, GC_section_lower, GC_section_upper, max_repeat_length):
    '''
    Accepts a list of three-letter residues or DNA sequence.
    Returns a list of resulting sequences with filtered GC-content and repetitive sequences.
    '''
    # Convert the DNA sequence string to protein list if necessary:
    if isinstance(seq, str):
        seq = DNA_to_protein(seq)
    # Initialise results list:
    results = []
    attempted_sequences = 0
    # While the number of results asked for is not reached, generate random sequences:
    while len(results) < no_of_results:
        potential_result = random_codon_assembly(seq)
        attempted_sequences += 1
        # Increments to notify of progress:
        increment = 1000
        if attempted_sequences % increment == 0:
            #if attempted_sequences != increment:
                #sys.stdout.flush()
            sys.stdout.write("\r" + str(attempted_sequences) + " sequences tested.")
            sys.stdout.flush()
        # Set a limit on number of random sequences generated:
        if attempted_sequences == max_rand_seq:
            sys.stdout.write("\r" + str(attempted_sequences) + " sequences tested. " + \
                             str(len(results)) + " sequences found.")
            break
        # Test if the random sequence has the desired GC content:
        if (GC_lower < GC_perc(potential_result) < GC_upper):
            # Test if each stretch of n bases is within the GC% threshold:
            if GC_sections(potential_result, GC_section_len, GC_section_lower, GC_section_upper) == True:
                # Test if there are no repetitive sequences longer than that specified:
                if DNA_repeats_filter(potential_result, max_repeat_length) == True:
                    results.append(potential_result)
    return(results)