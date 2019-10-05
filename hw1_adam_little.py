import re
from itertools import groupby

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

    loc = -1

    for a_acid in range(0, len(frame)):
        loc += 1
        site = frame[a_acid:a_acid+2]
        if site == 'KK' or site == 'KR' or site == 'RK' or site == 'RR':
            output_file.write("Double Basic site {} found at index: {}:{} \n".format(site, loc, loc+2))

def clevage_sites(frames):
    f1, f2, f3 = frames

    output_file.write("First reading frame:" + '\n')
    output_file.write(f1 + '\n')
    double_basic(f1)

    output_file.write("Second reading frame:" + '\n')
    output_file.write(f2 + '\n')
    double_basic(f2)

    output_file.write("Third reading frame:" + '\n')
    output_file.write(f3 + '\n')
    double_basic(f3)

output_file = open('translation_adam_little.fasta', 'w')

with open('pa1.fasta') as file:
    line = (x[1] for x in groupby(file, lambda l: l[0] == '>'))
    for header in line:
        header = next(header).strip()
        seq = ''.join(s.strip() for s in next(line))
        output_file.write(header + '\n')
        clevage_sites(dna_to_protein(seq))
