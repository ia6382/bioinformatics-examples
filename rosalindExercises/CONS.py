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
    f = open("40data.txt", "r")
    data = f.read()
    f.close

    S = parseFasta(data)

    symbols = {"A":0, "C":1, "G":2, "T":3}
    symbolsList = list(symbols.keys())

    profile = [[0]*len(S[0]) for i in range(4)]
    consesus = ""
    for i in range(len(S[0])):
        #build profile one position through all strings
        for s in S:
            profile[symbols[s[i]]][i] += 1

        #find most common nucleotide
        m = [0, 0]
        for j in range(4):
            if profile[j][i] > m[0]:
                m = [profile[j][i], j]
        consesus += symbolsList[m[1]]

    #write solution to file
    f = open("40result.txt", "w")
    f.write(consesus+"\n")
    for k, v in symbols.items():
        f.write(k+": ")
        for i in range(len(S[0])):
            f.write(str(profile[v][i])+" ")
        f.write("\n")
    f.close


