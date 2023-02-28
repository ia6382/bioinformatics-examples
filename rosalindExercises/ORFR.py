codonTable ={
"UUU":"F","CUU":"L","AUU":"I","GUU":"V",
"UUC":"F","CUC":"L","AUC":"I","GUC":"V",
"UUA":"L","CUA":"L","AUA":"I","GUA":"V",
"UUG":"L","CUG":"L","AUG":"M","GUG":"V",
"UCU":"S","CCU":"P","ACU":"T","GCU":"A",
"UCC":"S","CCC":"P","ACC":"T","GCC":"A",
"UCA":"S","CCA":"P","ACA":"T","GCA":"A",
"UCG":"S","CCG":"P","ACG":"T","GCG":"A",
"UAU":"Y","CAU":"H","AAU":"N","GAU":"D",
"UAC":"Y","CAC":"H","AAC":"N","GAC":"D",
"UAA":"*","CAA":"Q","AAA":"K","GAA":"E",
"UAG":"*","CAG":"Q","AAG":"K","GAG":"E",
"UGU":"C","CGU":"R","AGU":"S","GGU":"G",
"UGC":"C","CGC":"R","AGC":"S","GGC":"G",
"UGA":"*","CGA":"R","AGA":"R","GGA":"G",
"UGG":"W","CGG":"R","AGG":"R","GGG":"G"
}

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

def transcribeToRNA(dnaCoding):
    #RNA is a copy of dnaTemplate with U instead of T
    rna = dnaCoding.replace("T", "U")
    return rna

def translate(rna):
    protein = ""
    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]
        aminoAcid = codonTable[codon]
        protein += aminoAcid

    return protein


#read data
f = open("17data.txt", "r")
dna = f.read()
dnaComp = reverseComplement(dna)

#both strands of DNA can be used as a coding strand
rna = transcribeToRNA(dnaComp)
rnaComp = transcribeToRNA(dna)

#go through all possitions: 3 for different start nt *2 for complement
maxProtein = ""
for s in (rna, rnaComp):
    for frameStart in range(3):
        length = len(s) - (len(s)-frameStart)%3
        cutS = s[frameStart:length]

        #split tranlated s at stop codons
        acids = translate(cutS)
        potentialProteins = acids.split("*")
        
        #cut of acids before start acid M (=codon AUG)
        proteins = []
        for p in potentialProteins:
            for (j, a) in enumerate(p):
                if a == "M":
                    proteins.append(p[j:])

        #find max
        for p in proteins:
            if len(p) > len(maxProtein):
                maxProtein = p

#write solution to file
f = open("17result.txt", "w")
f.write(maxProtein)
f.close