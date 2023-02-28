def parse(data):
    lines = data.split("\n")

    return (lines)

def reverseComplement(dna):
    dnaC = ""
    for nt in dna[::-1]:
        if nt == "A":
            dnaC += "T"
        elif nt == "C":
            dnaC += "G"
        elif nt == "G":
            dnaC += "C"
        else:
            dnaC += "A"
    
    return dnaC


if __name__ == "__main__":
    #read file
    f = open("32data.txt", "r")
    data = f.read()
    f.close

    #parse input data into set
    dnaStrings = parse(data)
    k = len(dnaStrings[0])-1

    #create set S U Src
    S = set()
    [S.add(i) for i in dnaStrings]
    Src = set()
    [S.add(reverseComplement(i)) for i in dnaStrings]
    SUSrc = S | Src

    #connect nodes - create adjacency list
    al = []
    for i in SUSrc:
        node1 = i[0:k]
        node2 = i[1:k+1]
        al.append((node1, node2))

    #order list
    al.sort(key = lambda x: x[0])

    #write result
    f = open("32result.txt", "w")
    [f.write("("+i+", "+j+")\n") for (i,j) in al]
    f.close