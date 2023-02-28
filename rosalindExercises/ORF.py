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
"TAA":"s*","CAA":"Q","AAA":"K","GAA":"E",
"TAG":"s*","CAG":"Q","AAG":"K","GAG":"E",
"TGT":"C","CGT":"R","AGT":"S","GGT":"G",
"TGC":"C","CGC":"R","AGC":"S","GGC":"G",
"TGA":"s*","CGA":"R","AGA":"R","GGA":"G",
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


f = open("18data.txt", "r")
data = f.read()

#split FASTA format
fasta = data.split(">")[1:]
fasta = fasta[0]
label = fasta[:13]
dna = fasta[13:]
dna = dna.replace('\n', '')

#go through all possitions: 3 for different start nt *2 for complement
proteins = set() #all possible proteins
for s in (dna, reverseComplement(dna)):
    for frameStart in range(3):
        length = len(s) - (len(s)-frameStart)%3
        cutS = s[frameStart:length]

        #split tranlated s at stop codons (hack: stop codod marked as *s)
        acids = translate(cutS)
        potentialProteins = acids.split("*")
        
        #cut of acids before start acid M (=codon AUG)
        for p in potentialProteins:
            if p[-1] != "s": #to make sure it has a stop codon at the end
                continue
            for (j, a) in enumerate(p):
                if a == "M":
                    proteins.add(p[j:-1]) #ignore the last letter

#write solution to file
f = open("18result.txt", "w")
[f.write(p+"\n") for p in proteins]
f.close
'SHVANSGYMGMTPRLGLESLLE*A*MIRVASQ'