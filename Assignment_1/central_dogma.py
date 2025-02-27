import random
from tabulate import tabulate

MIN_PROTEIN_LENGTH = 25


def generate_DNA_sequence(length):
    """Generate a random DNA sequence of a given length"""
    sequence = ""
    nucleotides = ['A', 'T', 'C', 'G']
    for i in range(length):
        sequence += random.choice(nucleotides)
    return sequence

def check_nucleotide(sequence):
    """Check if a sequence is a valid DNA sequence"""
    valid_nucleotides = ['A', 'T', 'C', 'G']
    for nucleotide in sequence:
        if nucleotide not in valid_nucleotides:
            print(nucleotide)
            raise ValueError("Invalid nucleotide detected.")
    return True

def clean_sequence(sequence):
    '''Remove spaces and newlines from a DNA sequence'''
    sequence = sequence.upper()
    sequence = sequence.replace(" ", "")
    sequence = sequence.replace("\n", "")
    return sequence

def check_length(sequence):
    if len(sequence) % 3 != 0:
        print("Warning: The DNA sequence length is not divisible by 3.")
        # remove the last base
        sequence = sequence[:len(sequence) - (len(sequence) % 3)]
    return sequence
    
def complement(sequence):
    """Return the complement of a DNA sequence"""
    complement = ""
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    for nucleotide in sequence:
        if nucleotide in complement_dict:
            complement += complement_dict[nucleotide]
        else: 
            return "Invalid sequence"
    return complement

def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence"""
    return complement(sequence)[::-1]

def transcription(sequence):
    """Return the RNA transcript of a DNA sequence"""
    return sequence.replace('T','U')

def chunk_aa(sequence):
    """Return a list of codons from an RNA sequence"""
    codons = []
    for i in range(len(sequence)//3):
        codons.append(sequence[i*3:i*3+3])
    return codons


def translation(codons):
    """Return the protein sequence from a list of codons"""
    genetic_code = {
        "GCA": "Ala", "GCC": "Ala", "GCG": "Ala", "GCU": "Ala", #Alanine
        "AGA": "Arg", "AGG": "Arg", "CGA": "Arg", "CGC": "Arg", "CGG": "Arg", "CGU": "Arg", #Arginine
        "GAC": "Asp", "GAU": "Asp", #Aspartic acid
        "AAC": "Asn", "AAU": "Asn", #Asparagine
        "UGC": "Cys", "UGU": "Cys", #Cysteine
        "GAA": "Glu", "GAG": "Glu", #Glutamic Acid
        "CAA": "GIn", "CAG": "GIn", #Glutamine
        "GGA": "Gly", "GGC": "Gly", "GGU": "Gly", "GGG": "Gly", ##Glycine
        "CAC": "His", "CAU": "His", #Histidine
        "AUA": "IIe", "AUC": "IIe", "AUU": "IIe", #Isoleucine
        "UUA": "Leu", "UUG": "Leu", "CUA": "Leu", "CUC": "Leu", "CUG": "Leu", "CUU": "Leu", #Leucine
        "AAA": "Lys", "AAG": "Lys", #Lysine
        "AUG": "Met", #Methionine
        "UUC": "Phe", "UUU": "Phe", #Phenylalanine
        "CCA": "Pro", "CCC": "Pro", "CCG": "Pro", "CCU": "Pro", # Proline
        "AGC": "Ser", "AGU": "Ser", "UCA": "Ser", "UCC": "Ser", "UCG": "Ser", "UCU": "Ser", #Serine
        "ACA": "Thr", "ACC": "Thr", "ACG": "Thr", "ACU": "Thr", #Threonine
        "UGG": "Trp", #Tryptophan
        "UAC": "Tyr", "UAU": "Tyr", #Tyrosine
        "GUA": "Val", "GUC": "Val", "GUG": "Val", "GUU": "Val", #Valine
        "UAA": "Stop", "UAG": "Stop", "UGA": "Stop" #Termination codon
                  
    }

    protein = []
    for codon in codons:
        if len(codon) != 3:
            print(f"Debug: Invalid codon length: {codon}")
            return "Invalid sequence"
        
        if codon in genetic_code:
            if genetic_code[codon] == "Stop":
               # stop translation
                break
            protein.append(genetic_code[codon])
    
    return '-'.join(protein)


def dna_to_protein(sequence):
    """Return the protein sequence from a DNA sequence"""
    sequence = clean_sequence(sequence)
    
    # Validate nucleotide sequence
    try:
        check_nucleotide(sequence)
    except ValueError as e:
        print(f"Error: {e}")
        return None  # Return None if invalid

    sequence = check_length(sequence)  # This should still work as intended

    # Complement DNA
    complement_seq = complement(sequence)
    print(f"Complement DNA: {complement_seq}")

    # Transcription (Convert to mRNA)
    mrna = transcription(complement_seq)
    print(f"mRNA sequence: {mrna}")

    # Split mRNA into codons
    codons = chunk_aa(mrna)
    print(f"Codons: {codons}")

    # Translation to protein sequence
    aaseq = translation(codons)
    print(f"Protein sequence: {aaseq}")

    return aaseq


# 6-frame translation
def find_protein(sequence, strand, frame):
    """Find proteins in a given sequence strand and frame, filtering short proteins."""
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    proteins = {}
    
    seq_length = len(sequence)
    
    i = frame  # Start at the given frame

    while i < seq_length - 2:
        codon = sequence[i:i+3] # Loops through the sequence, reading 3 nucleotides at a time
        
        if codon == start_codon:  # Start translation
            start = i+1  # store the start index
            protein_seq = []
            
            for j in range(i, seq_length - 2, 3): # Iterate through the sequence, 3 nu at a time
                codon = sequence[j:j+3]
                if codon in stop_codons:  # Stop translation
                    end = j+3  # store the end index
                    length = (end - start) // 3 # Calculate protein length
                    if length >= MIN_PROTEIN_LENGTH:  # Filter out short proteins
                        if end not in proteins or length > proteins[end][2]:
                            proteins[end]  = (start, end, length, strand, frame+1, ''.join(protein_seq))
                    
                    break  # Move to the next start codon search
                protein_seq.append(codon)

        i += 3  # Move to next codon

    return list(proteins.values())


def six_frame_translation(sequence):
    """Perform 6-frame translation and detect proteins."""
    sequence = sequence.upper()
    
    frames = []
    
    # Forward frames
    for frame in range(3):
        frames.extend(find_protein(sequence, 'forward', frame))
    
    # Reverse complement
    reverse_seq = reverse_complement(sequence)
    seq_length = len(sequence)
    
    # Reverse frames
    for frame in range(3):
        rev_protein = find_protein(reverse_seq, 'reverse', frame)
        
        # Convert start/end to original coordinates
        for i in range(len(rev_protein)):
            start, end, length, strand, frame, protein = rev_protein[i]
            rev_protein[i] = (seq_length - end + 1, seq_length - start + 1, length, strand, frame, protein)
        
        frames.extend(rev_protein)
    
    return frames


def report_predicted_proteins(sequence):
    """Print predicted proteins in a table."""
    proteins = six_frame_translation(sequence)
    
    # Report proteins in table format
    if proteins:
        headers = ["Start", "End", "Length (aa)", "Strand", "Frame", "Sequence"]
        print(tabulate(proteins, headers=headers, tablefmt="pretty"))

        # Calculate longest and shortest proteins
        lengths = [info[2] for info in proteins]
        longest = max(lengths)
        shortest = min(lengths)

        print(f"\nLongest protein length: {longest} amino acids")
        for protein in proteins:
            if protein[2] == longest:
                print(f"Longest protein: Start={protein[0]}, End={protein[1]}, Length={protein[2]}, Strand={protein[3]}, Frame={protein[4]}")

        print(f"\nShortest protein length: {shortest} amino acids")
        for protein in proteins:
            if protein[2] == shortest:
                print(f"Shortest protein: Start={protein[0]}, End={protein[1]}, Length={protein[2]}, Strand={protein[3]}, Frame={protein[4]}")
    else:
        print("No valid protein found.")


    
    

    


