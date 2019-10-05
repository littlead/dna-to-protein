import re

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

    #Start codon
    'ATG': '?',

    # Stop codon
    'TAA': '!', 'TGA': '!', 'TAG': '!',
}

# General function used to tranlate each ORF
def protein(sequence, rf, start, stop):
    for base in range(start, stop, 3):
        codon = sequence[base:base+3]
        rf.append(amino[codon])

def dna_to_protein(sequence):
    seq = ''.join(sequence)
    # One amino acid sequence for each ORF
    rf1 = []
    rf2 = []
    rf3 = []
    if len(seq) % 3 == 0:
        protein(seq, rf1, 0, len(sequence))
        protein(seq, rf2, 1, len(sequence)-2)
        protein(seq, rf3, 2, len(sequence)-1)
    elif len(seq) % 3 == 1:
        protein(seq, rf1, 0, len(sequence)-1)
        protein(seq, rf2, 1, len(sequence))
        protein(seq, rf3, 2, len(sequence)-2)
    elif len(seq) %3 == 2:
        protein(seq, rf1, 0, len(sequence)-2)
        protein(seq, rf2, 1, len(sequence)-1)
        protein(seq, rf3, 2, len(sequence))
    else:
        print("Something weird is happening")

    frame1 = ''.join(rf1)
    frame2 = ''.join(rf2)
    frame3 = ''.join(rf3)

    return frame1, frame2, frame3

def double_basic(frame):

    loc = 0

    for a_acid in range(0, len(frame)):
        loc += 1
        site = frame[a_acid:a_acid+2]
        if site == 'KK' or site == 'KR' or site == 'RK' or site == 'RR':
            print("Double Basic site {} found at point: {}".format(site, loc))

def clevage_sites(frames):
    f1, f2, f3 = frames

    print("First reading frame:")
    print(f1)
    double_basic(f1)

    print("Second reading frame:")
    print(f2)
    double_basic(f2)

    print("Third reading frame:")
    print(f3)
    double_basic(f3)

sequence = []
seq = ''

with open('pa1.fasta', 'r') as file:

    lines = []

    for line in file:
        lines.append(line.rstrip())

    for line in lines:
        if line.startswith('>'):
            if sequence:
                seq = ''.join(sequence)
                frames = dna_to_protein(seq)
                clevage_sites(frames)
                sequence = []
                seq = ''
                print(line)
            else:
                print(line)
        else:
            sequence.append(line.rstrip())
