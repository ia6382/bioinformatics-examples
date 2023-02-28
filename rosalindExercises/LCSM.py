from collections import defaultdict

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

def LCS(s, t):
    lcsSet = set()
    lcsLen = 0

    #init matrix
    s = "-"+s
    t = "-"+t
    m = len(s)
    n = len(t)
    M = [[0]*n for i in range(m)]

    #fill matrix
    for i in range(1, m):
        for j in range(1, n):
            if s[i] == t[j]:
                M[i][j] = M[i-1][j-1] + 1
                if lcsLen < M[i][j]:
                    lcsLen = M[i][j]
                    lcs = backtrack(M, i, j, s)
                    lcsSet = {lcs}
                elif lcsLen == M[i][j]:
                    lcs = backtrack(M, i, j, s)
                    lcsSet.add(lcs)
    return lcsSet

def backtrack(M, i, j, s):
    lcs = ""
    while M[i][j] != 0:
        lcs += s[i]
        i -= 1
        j -= 1

    return lcs[::-1]

if __name__ == "__main__":
    #read file
    f = open("43data.txt", "r")
    data = f.read()
    f.close

    S = parseFasta(data)

    #longest common substring (LCS) of first two sequences
    maxLCSset = LCS(S[0], S[1])
    maxLCSlen = len(next(iter(maxLCSset))) #velikost prvega elementa

    #for all next squences
    for i in range(2, len(S)):
        currmaxLCSset = set()
        currmaxLCSlen = maxLCSlen

        #LCS' of all lcs' until now and next sequence
        for lcs in maxLCSset:
            lcsset = LCS(lcs, S[i])
            
            #remember only bigger or equaly long ones
            if len(lcsset) > 0:
                lcslen = len(next(iter(lcsset)))
                if lcslen > currmaxLCSlen:
                    currmaxLCSset = lcsset
                elif lcslen == currmaxLCSlen:
                    currmaxLCSset = currmaxLCSset | lcsset #union
        
        #updete LCS until this sequence
        maxLCSset = currmaxLCSset
        maxLCSlen = currmaxLCSlen

    print(maxLCSset)
    