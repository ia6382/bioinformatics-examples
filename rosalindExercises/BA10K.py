import math

def parse(data):
    lines = data.split("\n")

    #number of iterations
    numOfIter = int(lines[0])

    #emmisions (string of symbols)
    emmisions = lines[2]

    symbols = {}
    for i, s in enumerate(lines[4].split("\t")):
        symbols[s] = i
    m = len(symbols)

    #states
    states = {}
    for i, s in enumerate(lines[6].split("\t")):
        states[s] = i
    n = len(states)

    #probability state matrix
    pmatrixS = [[0]*n for i in range(n)]
    for i in range(n):
        line = lines[i+9].split("\t")
        for j in range(n):
            pmatrixS[i][j] = float(line[j+1])

    #probability emmisions matrix
    startIndexE = 11+len(states)
    pmatrixE = [[0]*m for i in range(n)]
    for i in range(n):
        line = lines[i+startIndexE].split("\t")
        for j in range(m):
            pmatrixE[i][j] = float(line[j+1])

    return (numOfIter, emmisions, symbols, states, pmatrixS, pmatrixE)

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

def gamma(f, b, k, i, n):
    g1 = f[i][k]*b[i][k]
    g2 = sum(f[i][l]*b[i][l] for l in range(n))
    
    return g1/g2

def eta(f, b, pmatrixS, pmatrixE, k, l, i, n, xi1):
    e1 = f[i][k] * pmatrixS[k][l] * b[i+1][l] * pmatrixE[l][xi1]
    sume = 0
    for y in range(n):
        for w in range(n):
            sume += f[i][y] * pmatrixS[y][w] * b[i+1][w] * pmatrixE[w][xi1]

    return e1/sume


if __name__ == "__main__":
    #read file
    f = open("31data.txt", "r")
    data = f.read()
    f.close

    #parse input data
    (numOfIter, emmisions, symbols, states, pmatrixS, pmatrixE) = parse(data)
    #symbol list for updating
    symbolList = sorted(list(symbols.items()), key=lambda x: x[1])
    symbolList = [x[0] for x in symbolList]
    #state list for writing result
    statesList = sorted(list(states.items()), key=lambda x: x[1])
    statesList = [x[0] for x in statesList]

    m = len(emmisions)
    n = len(states)

    #initial state distribution
    pi = [1/n for i in range(n)]

    for iteration in range(numOfIter):
        #forward, backward
        F = forward(emmisions, symbols, states, pmatrixS, pmatrixE, pi)
        B = backward(emmisions, symbols, states, pmatrixS, pmatrixE)

        #update initial state distribution
        #for i in range(n):
        #    pi[i] = gamma(F, B, i, 0, n)

        #update pmatrixS
        upMatS = [[0]*n for i in range(n)]
        for k in range(n):
            for l in range(n):
                a1 = sum(eta(F, B, pmatrixS, pmatrixE, k, l, i, n, symbols[emmisions[i+1]]) for i in range(m-1))
                a2 = sum(gamma(F, B, k, i, n) for i in range(m-1))

                upMatS[k][l] = a1/a2

        #update matrixE
        upMatE = [[0]*len(symbols) for i in range(n)]
        for k in range(n):
            for j in range(len(symbols)):
                e1 = sum((emmisions[i] == symbolList[j])*gamma(F, B, k, i, n) for i in range(m))
                e2 = sum(gamma(F, B, k, i, n) for i in range(m))

                upMatE[k][j] = e1/e2

        pmatrixS = upMatS
        pmatrixE = upMatE
    
    #write result
    f = open("31result.txt", "w")
    for i in range(n):
        f.write(statesList[i]+"\t")
    for i in range(n):
        f.write("\n"+statesList[i] +"\t")
        for j in range(n):
            f.write(str(round(pmatrixS[i][j], 3))+"\t")
    f.write("\n--------\n\t")
    for i in range(len(symbols)):
        f.write(symbolList[i]+"\t")
    for i in range(n):
        f.write("\n"+statesList[i]+"\t")
        for j in range(len(symbols)):
            f.write(str(round(pmatrixE[i][j], 3))+"\t")
    f.close

    """
    #forward-backward algorithm
    #combine into posterior probability matrix 
    pxF = sum(F[m-1])
    pxB = sum(pi[i] * pmatrixE[i][symbols[emmisions[0]]] * B[0][i] for i in range(0,n))

    P = [[0]*n for i in range(m)]
    for i in range(m):
        for j in range(n):
            P[i][j] = (F[i][j]*B[i][j])/pxF
    """