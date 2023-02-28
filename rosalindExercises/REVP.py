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
    f = open("41data.txt", "r")
    data = f.read()
    f.close

    dna = parseFasta(data)
    dna = dna[0]
    #dnaC = reverseComplement(dna)

    f = open("41result.txt", "w")
    for i in range(len(dna)):
        for j in range(4,13):
            if i+j <= len(dna):
                slice = dna[i:i+j]
                sliceC = reverseComplement(slice)

                if slice == sliceC:
                    f.write(str(i+1)+" "+str(j)+"\n")
    f.close