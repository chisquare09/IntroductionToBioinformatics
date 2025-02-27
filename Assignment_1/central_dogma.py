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
        "GCA": "A", "GCC": "A", "GCG": "A", "GCU": "A", #Alanine
        "AGA": "R", "AGG": "R", "CGA": "R", "CGC": "R", "CGG": "R", "CGU": "R", #Arginine
        "GAC": "D", "GAU": "D", #Aspartic acid
        "AAC": "N", "AAU": "N", #Asparagine
        "UGC": "C", "UGU": "C", #Cysteine
        "GAA": "E", "GAG": "E", #Glutamic Acid
        "CAA": "Q", "CAG": "Q", #Glutamine
        "GGA": "G", "GGC": "G", "GGU": "G", "GGG": "G", ##Glycine
        "CAC": "H", "CAU": "H", #Histidine
        "AUA": "I", "AUC": "I", "AUU": "I", #Isoleucine
        "UUA": "L", "UUG": "L", "CUA": "L", "CUC": "L", "CUG": "L", "CUU": "L", #Leucine
        "AAA": "K", "AAG": "K", #Lysine
        "AUG": "M", #Methionine
        "UUC": "F", "UUU": "F", #Phenylalanine
        "CCA": "P", "CCC": "P", "CCG": "P", "CCU": "P", # Proline
        "AGC": "S", "AGU": "S", "UCA": "S", "UCC": "S", "UCG": "S", "UCU": "S", #Serine
        "ACA": "T", "ACC": "T", "ACG": "T", "ACU": "T", #Threonine
        "UGG": "W", #Tryptophan
        "UAC": "Y", "UAU": "Y", #Tyrosine
        "GUA": "V", "GUC": "V", "GUG": "V", "GUU": "V", #Valine
        "UAA": "*", "UAG": "*", "UGA": "*" #Termination codon
                  
    }

    protein = []
    for codon in codons:
        if len(codon) != 3:
            print(f"Debug: Invalid codon length: {codon}")
            return "Invalid sequence"
        
        if codon in genetic_code:
            if genetic_code[codon] == "*":
               # stop translation
                break
            protein.append(genetic_code[codon])
    
    return ''.join(protein)


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
    proteins = []
    
    seq_length = len(sequence)
    i = frame  # Start at the given frame

    while i < seq_length - 2:
        codon = sequence[i:i+3]
        
        if codon == start_codon:  # Found a start codon
            start = i + 1  # 1-based index
            protein_seq = []
            
            for j in range(i, seq_length - 2, 3):
                codon = sequence[j:j+3]
                if codon in stop_codons:  # Stop translation
                    end = j + 3  # 1-based inclusive index
                    length = (end - start) // 3
                    
                    if length >= MIN_PROTEIN_LENGTH:  # Filter out short proteins
                        proteins.append((start, end, length, strand, frame + 1, "".join(protein_seq)))
                    
                    break  # Move to the next start codon search
                protein_seq.append(codon)
        
        i += 3  # Move to next codon

    return proteins


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


    
    

    


