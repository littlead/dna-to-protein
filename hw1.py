
data = 'ATAAGC'

amino = {
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'TTA': 'L', 'TTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TTT': 'F', 'TTC': 'F',
    'TGG': 'W',
    'TAT': 'Y', 'TAC': 'Y',
    'GAT': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    'CAT': 'H', 'CAC': 'H',
    'AAA': 'K', 'AAG': 'K',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'TGT': 'C', 'TGC': 'C',
    'ATG': 'M',
    'AAT': 'N', 'AAC': 'N',
    'CAA': 'Q', 'CAG': 'Q',
    'ATG': 'START',
    'TAA': 'STOP', 'TGA': 'STOP', 'TAG': 'STOP',
}

# General function used to tranlate each ORF
def protein(sequence, orf, start, stop):
    for base in range(start, stop, 3):
        codon = sequence[base:base+3]
        orf.append(amino[codon])

def dna_to_protein(sequence):
    
    # One amino acid sequence for each ORF
    orf1 = []
    orf2 = []
    orf3 = []

    protein(sequence, orf1, 0, len(sequence))
    protein(sequence, orf2, 1, len(sequence)-2)
    protein(sequence, orf3, 2, len(sequence)-1)

    # Print to console for testing
    print(orf1, orf2, orf3)

    return orf1, orf2, orf3

# Call Function
dna_to_protein(data)
