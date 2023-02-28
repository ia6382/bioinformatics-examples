dnaCodonTable ={
"TTT":"F","CTT":"L","ATT":"I","GTT":"V",
"TTC":"F","CTC":"L","ATC":"I","GTC":"V",
"TTA":"L","CTA":"L","ATA":"I","GTA":"V",
"TTG":"L","CTG":"L","ATG":"M","GTG":"V",
"TCT":"S","CCT":"P","ACT":"T","GCT":"A",
"TCC":"S","CCC":"P","ACC":"T","GCC":"A",
"TCA":"S","CCA":"P","ACA":"T","GCA":"A",
"TCG":"S","CCG":"P","ACG":"T","GCG":"A",
"TAT":"Y","CAT":"H","AAT":"N","GAT":"D",
"TAC":"Y","CAC":"H","AAC":"N","GAC":"D",
"TAA":"s","CAA":"Q","AAA":"K","GAA":"E",
"TAG":"s","CAG":"Q","AAG":"K","GAG":"E",
"TGT":"C","CGT":"R","AGT":"S","GGT":"G",
"TGC":"C","CGC":"R","AGC":"S","GGC":"G",
"TGA":"s","CGA":"R","AGA":"R","GGA":"G",
"TGG":"W","CGG":"R","AGG":"R","GGG":"G"
}

def translate(dna):
    protein = ""
    for i in range(0, len(dna), 3):
        codon = dna[i:i+3]
        aminoAcid = dnaCodonTable[codon]
        protein += aminoAcid

    return protein

def reverseComplement(dna):
    dnaC = ""
    for nt in dna[::-1]:
        if nt == "A":
            dnaC += "T"
        elif nt == "C":
            dnaC += "G"
        elif nt == "G":
            dnaC += "C"
        else:
            dnaC += "A"
    
    return dnaC

def parseFasta(data):
    #split FASTA format
    fasta = data.split(">")[1:]
    S = []
    for i in fasta:
        label = i[:14]
        strand = i[14:]
        strand = strand.replace('\n', '')
        S.append(strand)

    return S


if __name__ == "__main__":
    #read file
    f = open("42data.txt", "r")
    data = f.read()
    f.close

    S = parseFasta(data)
    dna = S[0]
    introns = S[1:]

    for i in introns:
        dna = dna.replace(i, "")

    #dnaC = reverseComplement(dna)
    protein = translate(dna)[:-1]
    print(protein)