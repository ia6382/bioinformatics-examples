def composition(sequence, k):
    kmers = []
    l = len(sequence)
    for i in range(l):
        if i+k <= l:
            kmers.append(sequence[i:i+k])
        else:
            #if we need to loop around create kmer character by character
            kmer = ""
            for j in range(k):
                kmer += sequence[(i+j)%l]
            kmers.append(kmer)

    return kmers

def deBrujinGraph(dnaStrings):
    k = len(dnaStrings[0])-1
    SUSrc = [] #na rosalindu mora biti set, samo potem ne handla repeats
    [SUSrc.append(i) for i in dnaStrings]
    
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
    candidates = []
    stuck = True
    for i, nextConn in enumerate(al):
        if node == nextConn[0]: #we found an outgoing connection
            stuck = False
            candidates.append(i)

    return stuck, candidates


if __name__ == "__main__":
    #read file
    f = open("34data.txt", "r")
    data = f.read()
    dnaStrings = data.split("\n")
    f.close

    #parse input data into set
    #dnaStrings = composition(data,2)

    #adjecency list of de Brujin graph
    al = deBrujinGraph(dnaStrings)
    
    #Hierholzer algorithm for Euler cycle with queue
    #find a cycle for every branch that we come across
    Q = []
    A = []
    C = []

    conn = al[0]
    queue = [conn[0]]
    al.pop(0)

    Q.append(queue)
    A.append(al)
    C.append(conn)

    iterations = len(al)
    while(iterations > 0):
        for j in range(len(Q)):
            stuck, candidates = findNextConn(C[j][1], A[j])
            while stuck:
                lastEl = Q[j].pop()
                Q[j] = [lastEl]+Q[j]
                lastEl = Q[j][-1]
                stuck, candidates = findNextConn(lastEl, A[j])

            origQueue = Q[j].copy() #remeber data before changing it 
            origAl = A[j].copy()
            for e, i in enumerate(candidates):
                if e > 0: #check all branches
                    newQueue = origQueue.copy()
                    newAl = origAl.copy()

                    newConn = newAl[i]
                    newQueue = [newConn[0]]+newQueue
                    newAl.pop(i)

                    Q.append(newQueue)
                    A.append(newAl)
                    C.append(newConn)
                else:
                    newConn = A[j][i]
                    Q[j] = [newConn[0]]+Q[j]
                    A[j].pop(i)
                    C[j] = newConn


        iterations -= 1

    cycles = []
    for i in Q:
        #reverse queues
        queue = i
        queue = queue[::-1]

        #get cyclic genome
        genome = ""
        for i in queue:
            genome += i[-1]

        #add only unique cycles
        cycles.append(genome)
    
    #keep only unique cycles
    for i in range(len(cycles)-1, -1, -1):
        for j in range(i-1, -1, -1):
            c1 = cycles[i]
            c2 = cycles[j]
            C = c2+c2
            if c1 in C:
                cycles.pop(i)
                break
    
    #v rosalindu bi moral se zaciklati vsakega da se zacne z prvim kmerom

    [print(i) for i in cycles]
    print(len(cycles))