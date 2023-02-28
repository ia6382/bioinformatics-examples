import math

def parse(data):
    lines = data.split("\n")

    #emmisions (string of symbols)
    emmisions = lines[0]

    symbols = {}
    for i, s in enumerate(lines[2].split("\t")):
        symbols[s] = i
    m = len(symbols)

    #path (string of states)
    path = lines[4]

    states = {}
    for i, s in enumerate(lines[6].split("\t")):
        states[s] = i
    n = len(states)

    #infer probability state matrix
    pmatrixS = [[0]*n for i in range(n)]
    #count all types of transitions
    for i in range(len(path)-1):
        k = path[i]
        l = path[i+1]
        pmatrixS[states[k]][states[l]] += 1
    #divide by sum of row for that type => prob
    for i in range(n):
        rowsum = sum(pmatrixS[i])
        for j in range(n):
            pmatrixS[i][j] = pmatrixS[i][j]/ rowsum


    #probability emmisions matrix
    pmatrixE = [[0]*m for i in range(n)]
    for i in range(n):
        line = lines[i+9].split("\t")
        for j in range(m):
            pmatrixE[i][j] = float(line[j+1])

    return (emmisions, symbols, path, states, pmatrixS, pmatrixE)

def probPath(path, states, pmatrixS):
    #probability of path
    pPath = math.log(1)
    for i in range(len(path)-1):
        k = path[i]
        l = path[i+1]
        akl = pmatrixS[states[k]][states[l]]
        pPath += math.log(akl)

    pPath = math.exp(pPath)
    return(pPath)

def jointProbEmmisionsPath(emmisions, symbols, path, states, pmatrixS, pmatrixE):
    prob = math.log(1)
    for i in range(len(path)-1):
        k = path[i]
        l = path[i+1]
        akl = pmatrixS[states[k]][states[l]]
        
        x = emmisions[i]
        ekx = pmatrixE[states[k]][symbols[x]]

        prob += math.log(ekx*akl)


    k = path[len(path)-1]
    x = emmisions[len(path)-1]
    ekx = pmatrixE[states[k]][symbols[x]]

    prob += math.log(ekx)

    prob = math.exp(prob)
    return prob

def jointProbEmmisionsPath2(emmisions, symbols, path, states, pmatrixS, pmatrixE):
    prob = math.log(1)
    for i in range(len(path)):
        k = path[i]
        
        x = emmisions[i]
        ekx = pmatrixE[states[k]][symbols[x]]

        prob += math.log(ekx)
    
    return math.exp(prob)

if __name__ == "__main__":
    #read file
    f = open("28data.txt", "r")
    data = f.read()
    f.close

    (emmisions, symbols, path, states, pmatrixS, pmatrixE) = parse(data)

    pPath = probPath(path, states, pmatrixS)
    pEmmisionsPath = jointProbEmmisionsPath(emmisions, symbols, path, states, pmatrixS, pmatrixE)

    pEmmisionsCondPath = pEmmisionsPath / pPath
    print(pEmmisionsCondPath)
    
    #poenostavljena formula -> verjetnosti a se zaradi deljenjja pokrajsajo
    pSimplified = jointProbEmmisionsPath2(emmisions, symbols, path, states, pmatrixS, pmatrixE)
    print(pSimplified)
