
data = 'ATAAGC'

def dna_to_protein(sequence):

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


    orf1 = []
    orf2 = []
    orf3 = []

    for base in range(0, len(sequence), 3):
        codon = sequence[base:base+3]
        orf1.append(amino[codon])

    for base in range(1, len(sequence)-2, 3):
        codon = sequence[base:base+3]
        orf2.append(amino[codon])

    for base in range(2, len(sequence)-1, 3):
        codon = sequence[base:base+3]
        orf3.append(amino[codon])

    # Print to console for testing
    print(orf1, orf2, orf3)
    return orf1, orf2, orf3

# Call Function
dna_to_protein(data)
