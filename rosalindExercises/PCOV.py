def parse(data):
    lines = data.split("\n")

    return (lines)

def deBrujinGraph(dnaStrings):
    k = len(dnaStrings[0])-1
    SUSrc = set() #na rosalindu mora biti set, samo potem ne handla repeats
    [SUSrc.add(i) for i in dnaStrings]
    
    #connect nodes - create adjacency list
    al = []
    for i in SUSrc:
        node1 = i[0:k]
        node2 = i[1:k+1]
        al.append((node1, node2))

    #order list
    al.sort(key = lambda x: x[0])

    return al

def findNextConn(node, al):
    stuck = True
    for i, nextConn in enumerate(al):
        if node == nextConn[0]:
            stuck = False
            break

    return stuck, i


if __name__ == "__main__":
    #read file
    f = open("33data.txt", "r")
    data = f.read()
    f.close

    #parse input data into set
    dnaStrings = parse(data)
    
    #adjecency list of de Brujin graph
    al = deBrujinGraph(dnaStrings)

    #Hierholzer algorithm for Euler cycle with queue
    conn = al[0]
    queue = [conn[0]]
    al.pop(0)

    while(len(al) > 0):
        stuck, i = findNextConn(conn[1], al)
        while stuck:
            lastEl = queue.pop()
            queue = [lastEl]+queue
            lastEl = queue[-1]
            stuck, i = findNextConn(lastEl, al)
            
        conn = al[i]
        queue = [conn[0]]+queue
        al.pop(i)

    queue = queue[::-1]

    #get cyclic genome
    genome = ""
    for i in queue:
        genome += i[-1] 
    print(genome)

    """
    test v zvezku:
    01
    12
    34
    20
    13
    41
    """