import math

def parse(data):
    lines = data.split("\n")

    #emmisions (string of symbols)
    emmisions = lines[0]

    symbols = {}
    for i, s in enumerate(lines[2].split("\t")):
        symbols[s] = i
    m = len(symbols)

    path = lines[4]

    #states
    states = {}
    for i, s in enumerate(lines[6].split("\t")):
        states[s] = i
    n = len(states)

    return (emmisions, symbols, path, states)

def estimatePmatrixS(path, states):
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
        if rowsum > 0:
            for j in range(n):
                pmatrixS[i][j] = pmatrixS[i][j]/ rowsum

    return pmatrixS

def estimatePmatrixE(emmisions, symbols, path, states):
    m = len(states)
    n = len(symbols)

    #infer probability emssion state matrix
    pmatrixE = [[0]*n for i in range(m)]

    #count number of symbols x given state k
    for i in range(len(emmisions)):
        k = path[i]
        x = emmisions[i]

        pmatrixE[states[k]][symbols[x]] += 1
    
    #divide by sum of row for that type => prob
    for i in range(m):
        rowsum = sum(pmatrixE[i])
        if rowsum > 0:
            for j in range(n):
                pmatrixE[i][j] = pmatrixE[i][j]/ rowsum

    return pmatrixE


if __name__ == "__main__":
    #read file
    f = open("37data.txt", "r")
    data = f.read()
    f.close

    #parse input data
    (emmisions, symbols, path, states) = parse(data)

    m = len(emmisions)
    n = len(states)

    #symbol list for writing results
    symbolList = sorted(list(symbols.items()), key=lambda x: x[1])
    symbolList = [x[0] for x in symbolList]
    #state list for writing result
    statesList = sorted(list(states.items()), key=lambda x: x[1])
    statesList = [x[0] for x in statesList]

    #initial state distribution
    pi = [1/n for i in range(n)]

    #update parameters
    pmatrixS = estimatePmatrixS(path, states)
    pmatrixE = estimatePmatrixE(emmisions, symbols, path, states)

    #write result
    f = open("37result.txt", "w")
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
