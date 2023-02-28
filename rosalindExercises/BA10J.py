import math

def parse(data):
    lines = data.split("\n")

    #emmisions (string of symbols)
    emmisions = lines[0]

    symbols = {}
    for i, s in enumerate(lines[2].split("\t")):
        symbols[s] = i
    m = len(symbols)

    #states
    states = {}
    for i, s in enumerate(lines[4].split("\t")):
        states[s] = i
    n = len(states)

    #probability state matrix
    pmatrixS = [[0]*n for i in range(n)]
    for i in range(n):
        line = lines[i+7].split("\t")
        for j in range(n):
            pmatrixS[i][j] = float(line[j+1])

    #probability emmisions matrix
    startIndexE = 9+len(states)
    pmatrixE = [[0]*m for i in range(n)]
    for i in range(n):
        line = lines[i+startIndexE].split("\t")
        for j in range(m):
            pmatrixE[i][j] = float(line[j+1])

    return (emmisions, symbols, states, pmatrixS, pmatrixE)

def forward(emmisions, symbols, states, pmatrixS, pmatrixE, pi):
    #dynamic programming
    m = len(emmisions)
    n = len(states)
    F = [[0]*n for i in range(m)]

    #initialize: sum(F[prejsni]) = 1 (only begin state)
    for i in range(n):
        x = emmisions[0]
        eix = pmatrixE[i][symbols[x]]
        akj = pi[i] #1/len(states)
        F[0][i] = eix * 1 * akj 

    #traverse matrix
    for i in range(1, m):
        x = emmisions[i]
        for j in range(0, n):
            #from all states
            fs = []
            for k in range(0, n):
                akj = pmatrixS[k][j]
                fs.append(F[i-1][k] * akj)

            ejx = pmatrixE[j][symbols[x]]
            sumfs = sum(fs)
            F[i][j] = ejx * sumfs
    
    return F

def backward(emmisions, symbols, states, pmatrixS, pmatrixE):
    #dynamic programming
    m = len(emmisions)
    n = len(states)
    B = [[0]*n for i in range(m)]

    #initialize last row is 1
    for i in range(n):
        B[m-1][i] = 1

    #traverse matrix backwards
    for i in range(m-2, -1, -1):
        x = emmisions[i+1] #po simbolih gremo od [m-1 do 1], ker moramo nekako priditi v simbol 0
        for j in range(0, n):
            #for all states
            bs = []
            for k in range(0, n):
                ajk = pmatrixS[j][k]
                ekx = pmatrixE[k][symbols[x]]
                bs.append(ajk*ekx*B[i+1][k])

            sumbs = sum(bs)
            B[i][j] = sumbs
    
    return B


if __name__ == "__main__":
    #read file
    f = open("38data.txt", "r")
    data = f.read()
    f.close

    #parse input data
    (emmisions, symbols, states, pmatrixS, pmatrixE) = parse(data)

    m = len(emmisions)
    n = len(states)

    #state list for writing result
    statesList = sorted(list(states.items()), key=lambda x: x[1])
    statesList = [x[0] for x in statesList]

    #forward-backward algorithm
    #initial state distribution
    pi = [1/n for i in range(n)]

    #forward, backward
    F = forward(emmisions, symbols, states, pmatrixS, pmatrixE, pi)
    B = backward(emmisions, symbols, states, pmatrixS, pmatrixE)

    #combine into posterior probability matrix 
    pxF = sum(F[m-1])
    #pxB = sum(pi[i] * pmatrixE[i][symbols[emmisions[0]]] * B[0][i] for i in range(0,n))

    P = [[0]*n for i in range(m)]
    for i in range(m):
        for j in range(n):
            P[i][j] = (F[i][j]*B[i][j])/pxF

    #write result
    f = open("38result.txt", "w")
    for i in range(n):
        f.write(statesList[i]+"\t")
    for i in range(m):
        f.write("\n")
        for j in range(n):
            f.write(str(round(P[i][j], 4))+"\t")
    f.close