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

def pDistance(s1, s2):
    hammingDistance = 0
    for (i, j) in zip(s1, s2):
        if i != j:
            hammingDistance += 1
    
    return hammingDistance/len(s1) #ratio distance

if __name__ == "__main__":
    #read file
    f = open("19data.txt", "r")
    data = f.read()
    f.close

    S = parseFasta(data)

    #build distance matrix
    n = len(S)
    D = [[0]*n for i in range(n)]
    for i in range(n):
        for j in range(i):
            d = pDistance(S[i], S[j])
            D[i][j] = d
            D[j][i] = d

    #write solution to file
    f = open("19result.txt", "w")
    for i in range(n):
        for j in range(n):
            f.write(str(D[i][j]))
            if j < n-1:
                f.write(" ")
        if i < n-1:
            f.write("\n")
    f.close
            