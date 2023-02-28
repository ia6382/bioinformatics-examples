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

def editScore(bp1, bp2):
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
    f = open("22data.txt", "r")
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
        M[i][0] = (i,"down")

    for i in range(n):
        M[0][i] = (i, "right")

    for i in range(1, m):
        for j in range(1, n):
            m1 = (M[i-1][j][0] + editScore(s1[i], "-"), "down")
            m2 = (M[i][j-1][0] + editScore("-", s2[j]), "right")
            m3 = (M[i-1][j-1][0] + editScore(s1[i], s2[j]), "diag")
            M[i][j] = min(m1, m2, m3)

    #backtrack 
    i = m-1
    j = n-1
    optAls1 = ""
    optAls2 = ""
    while i != 0 and j != 0:
        direction = M[i][j][1]
        if direction == "down":
            optAls1 += s1[i]
            optAls2 += "-"
            i -= 1
        elif direction == "right":
            optAls1 += "-"
            optAls2 += s2[j]
            j -= 1
        else: #direction == "diag"
            optAls1 += s1[i]
            optAls2 += s2[j]
            i -= 1
            j -= 1

    optAls1 = optAls1[::-1]
    optAls2 = optAls2[::-1]
    #write solution to file
    f = open("22result.txt", "w")
    f.write(str(M[m-1][n-1][0])+"\n"+optAls1+"\n"+optAls2)
    f.close()