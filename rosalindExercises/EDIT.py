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

def editOperationCost(bp1, bp2):
    d = 0
    if bp1 == "-" or bp2 == "-":
        d = 1
    elif bp1 != bp2: #and bp1 != "-" and bp2 != "-":
        d = 1
    else: #bp1 == bp2
        d = 0

    return d

if __name__ == "__main__":
    #read file
    f = open("21data.txt", "r")
    data = f.read()
    f.close

    (s1, s2) = parseFasta(data)

    #dynamic programming global alingment cost
    s1 = "-"+s1
    s2 = "-"+s2
    m = len(s1)
    n = len(s2)
    M = [[0]*n for i in range(m)]

    for i in range(m):
        M[i][0] = i

    for i in range(n):
        M[0][i] = i

    for i in range(1, m):
        for j in range(1, n):
            m1 = M[i-1][j] + editOperationCost(s1[i], "-") #down
            m2 = M[i][j-1] + editOperationCost("-", s2[j]) #right
            m3 = M[i-1][j-1] + editOperationCost(s1[i], s2[j]) #diag
            M[i][j] = min(m1, m2, m3)

    print(M[m-1][n-1])