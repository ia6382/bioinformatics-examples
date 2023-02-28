def parse(data):
    genomes = []
    lines = data.split("\n")
    for i in range(0, len(lines)-1, 2):
        genomes.append((lines[i][2:], lines[i+1])) #(name, sequence)

    return (genomes)

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
    S = [] #na rosalindu mora biti set, samo potem ne handla repeats
    [S.append(i) for i in dnaStrings]
    
    #connect nodes - create adjacency list
    al = []
    for i in S:
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

def assemble(al):
    #Hierholzer algorithm for Euler cycle with queue
    # + keep track and find a cycle for every branch that we come across
    Q = [] #list of all queues we are keeping track of
    A = [] #list of all adjeceny list
    C = [] #list of all next connections

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

    #extract solutions (cycles)
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
    
    return cycles


if __name__ == "__main__":
    #read file
    f = open("genomes.txt", "r")
    data = f.read()
    f.close

    #parse input data into set
    genomes = parse(data)
    
    for (name, sequence) in genomes:
        print("\n"+name)

        #first roughly estimate required kmer length
        estK = 2
        for k in range(2,100):
            #split sequence into kmers of lenght k
            kmers = composition(sequence, k)

            #adjecency list of de Brujin graph
            al = deBrujinGraph(kmers)

            #branching in graph = multiple possible cycles = no unique solution in our test data
            branching = False
            nodes = set()
            for (fromNode, toNode) in al:
                if fromNode in nodes:
                    branching = True
                    break
                else:
                    nodes.add(fromNode)
            
            if branching == False:
                estK = k
                print("\t***estimated lenght of kmer: "+str(estK))
                break
        
        #shorten required kmer length by precisely finding all cycles
        assembledSequence = ""
        for k in range(estK, 2, -1):
            #split sequence into kmers of lenght k
            kmers = composition(sequence, k)

            #adjecency list of de Brujin graph
            al = deBrujinGraph(kmers)

            #find all possible cycles in graph
            cycles = assemble(al)

            #if there is more than one cycle it means solution is not unique
            if len(cycles) > 1:
                #check result if unique solution was found previous iteration              
                a = assembledSequence + assembledSequence #our string may be just shifted original string
                if sequence in a:
                    print("\tFound minimum lenght of kmers for unique reconstruction: "+str(k+1))
                    #print(assembledSequence)
                else:
                    print("\tOriginal sequence does not match assembled sequence."+str(k+1))
                    #print(assembledSequence)
                break
            else:
                assembledSequence = cycles[0]
    