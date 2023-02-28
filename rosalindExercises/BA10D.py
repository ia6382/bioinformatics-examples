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


if __name__ == "__main__":
    #read file
    f = open("36data.txt", "r")
    data = f.read()
    f.close

    #parse input data
    (emmisions, symbols, states, pmatrixS, pmatrixE) = parse(data)

    m = len(emmisions)
    n = len(states)

    #initial state distribution
    pi = [1/n for i in range(n)]

    #forward algorithm
    F = forward(emmisions, symbols, states, pmatrixS, pmatrixE, pi)

    #probability for emmited sequence given the HMM
    pxF = sum(F[m-1])
    print(pxF)