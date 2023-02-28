import math

def parse(data):
    lines = data.split("\n")
    
    #path - (string of states)
    path = lines[0]
    
    #states
    states = {}
    for i, s in enumerate(lines[2].split("\t")):
        states[s] = i
    n = len(states)

    #probability matrix
    pmatrix = [[0]*n for i in range(n)]
    for i in range(n):
        line = lines[i+5].split("\t")
        for j in range(n):
            pmatrix[i][j] = float(line[j+1])

    return (path, states, pmatrix)

if __name__ == "__main__":
    #read file
    f = open("27data.txt", "r")
    data = f.read()
    f.close

    (path, states, pmatrix) = parse(data)

    #probability of path
    pPath = math.log(1/len(states)) #starting probabilities are equal
    for i in range(len(path)-1):
        k = path[i]
        l = path[i+1]
        akl = pmatrix[states[k]][states[l]]
        pPath += math.log(akl)

    pPath = math.exp(pPath)
    print(pPath)

