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
    pmatrixE = [[0]*m for i in range(n)]
    for i in range(n):
        line = lines[i+11].split("\t")
        for j in range(m):
            pmatrixE[i][j] = float(line[j+1])

    return (emmisions, symbols, states, pmatrixS, pmatrixE)

if __name__ == "__main__":
    #read file
    f = open("29data.txt", "r")
    data = f.read()
    f.close

    #parsa input data
    (emmisions, symbols, states, pmatrixS, pmatrixE) = parse(data)
    #create list of states for backtracking
    statesList = sorted(list(states.items()), key=lambda x: x[1])
    statesList = [x[0] for x in statesList]

    #dynamic programming
    m = len(emmisions)
    n = len(states)
    V = [[0]*n for i in range(m)]

    #initialize: max(V[prejsni]) = 1 (from begin state)
    for i in range(n):
        x = emmisions[0]
        eix = pmatrixE[i][symbols[x]]
        akj = 1/len(states)
        V[0][i] = (eix * 1 * akj, "") 

    #traverse matrix
    for i in range(1, m):
        x = emmisions[i]
        for j in range(0, n):
            #from all states
            vs = []
            for k in range(0, n):
                akj = pmatrixS[k][j]
                vs.append((V[i-1][k][0] * akj, statesList[k]))

            ejx = pmatrixE[j][symbols[x]]
            maxvs = max(vs)
            V[i][j] = (ejx * maxvs[0], maxvs[1])

    #backtrack from max el in last row
    maxi = V[m-1].index(max(V[m-1]))
    maxel = V[m-1][maxi]
    path = statesList[maxi]

    i = m-1
    while i > 0:
        prevState = maxel[1]
        path += prevState

        i -= 1
        maxel = V[i][states[prevState]]
    
    path = path[::-1]
    print(path)
    #BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB