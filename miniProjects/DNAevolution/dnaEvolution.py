import numpy as np
import math
import time

from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.cluster import AgglomerativeClustering

def parseGenes(data):
    data = data.split("\n")
    genes = {}
    l = len(data)
    for i in range(0,l-8, 8):
        gene = data[i]
        s = data[i+1:i+8]
        species = []
        for j in s:
            j = j.split("\t")
            name = j[1]
            sequence = j[2]
            species.append((name, sequence))
        genes[gene] = species

    return genes

def editScore(bp1, bp2):
    d = 0
    if bp1 == "-" or bp2 == "-":
        d = -2
    elif bp1 != bp2: #and bp1 != "-" and bp2 != "-":
        d = -1
    else: #bp1 == bp2
        d = 1

    return d

def globalAlingment(s1, s2):
    #dynamic programming global alingment cost
    s1 = "-"+s1
    s2 = "-"+s2
    m = len(s1)
    n = len(s2)
    M = [[0]*n for i in range(m)]

    for i in range(m):
        M[i][0] = (i,"down")

    for i in range(n):
        M[0][i] = (i, "right")

    for i in range(1, m):
        for j in range(1, n):
            m1 = (M[i-1][j][0] + editScore(s1[i], "-"), "down")
            m2 = (M[i][j-1][0] + editScore("-", s2[j]), "right")
            m3 = (M[i-1][j-1][0] + editScore(s1[i], s2[j]), "diag")
            M[i][j] = max(m1, m2, m3)

    #backtrack 
    i = m-1
    j = n-1
    optAls1 = ""
    optAls2 = ""
    while i != 0 and j != 0:
        direction = M[i][j][1]
        if direction == "down":
            optAls1 += s1[i]
            optAls2 += "-"
            i -= 1
        elif direction == "right":
            optAls1 += "-"
            optAls2 += s2[j]
            j -= 1
        else: #direction == "diag"
            optAls1 += s1[i]
            optAls2 += s2[j]
            i -= 1
            j -= 1

    optAls1 = optAls1[::-1]
    optAls2 = optAls2[::-1]

    return (optAls1, optAls2)

def jukesCantor(s1, s2):
    #proportion of differences
    diff = 0
    for n1, n2 in zip(s1, s2):
        if n1 != n2: #or n1 == "-" or n2 == "-":
            diff += 1
    p = diff/len(s1)

    p = min(p, (3/4)) #limit value of p because of ln

    #JC formula
    eDistance = (-3/4) * math.log(1 - (3/4) * p)

    return eDistance

def distMatrix(genes, gene):
    species = genes[gene]
    l = len(species)

    D = np.zeros((l,l))
    for i in range(l):
        for j in range(i+1, l):
            s1 = species[i][1]
            s2 = species[j][1]

            (optAls1, optAls2) = globalAlingment(s1, s2)
            eDistance = jukesCantor(optAls1, optAls2)

            D[i][j] = eDistance
            D[j][i] = eDistance

    return D

def plot_dendrogram(model):
    """
    https://scikit-learn.org/stable/auto_examples/cluster/plot_agglomerative_dendrogram.html?highlight=hierarchical
    """
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack([model.children_, model.distances_,
                                      counts]).astype(float)

    # Plot the corresponding dendrogram
    labels = ['Bos taurus', 'Mus musculus', 'Homo sapiens', 'Danio rerio', 'Canis lupus familiaris', 'Rattus norvegicus', 'Macaca mulatta']
    dendrogram(linkage_matrix, orientation='left', labels=labels)

def phylogeneticTree(genes, gene):
    #distance matrix
    D = distMatrix(genes, gene)

    #phylogenetic tree UPGMA
    model = AgglomerativeClustering(affinity='precomputed', linkage="average", distance_threshold=0, n_clusters=None)
    model = model.fit(D)
    
    #draw dendrogram
    plot_dendrogram(model)

    plt.title("Gene "+gene+" Phylogenetic Tree")
    plt.ylabel("Animal Species")
    plt.xlabel("Evolutionary Distance")
    ax = plt.gca() #get current axis
    ax.tick_params(axis='x', which='major', labelsize=10)
    ax.tick_params(axis='y', which='major', labelsize=10)
    plt.tight_layout()
    plt.show()

def rankGenes(genes):
    """
    using average distance between groups
    """
    ranking = []

    for gene in genes.keys():
        species = genes[gene]
        s1 = species[3][1] #fish 

        #average distance between groups
        d = 0
        for i in [0, 1, 2, 4, 5, 6]: #mammals
            s2 = species[i][1]
            (optAls1, optAls2) = globalAlingment(s1, s2)
            eDistance = jukesCantor(optAls1, optAls2)
            d += eDistance
        d = d/6

        ranking.append((gene, d))
    
    #sort decreasing
    ranking = sorted(ranking, key=lambda x: x[1], reverse=True)
    return ranking

def rankGenes2(genes):
    """
    using heuristic (cut of ends so the strings are same lenght) consensus of mammals
    """
    ranking = []
    symbols = {"A":0, "C":1, "G":2, "T":3, "N":4}
    symbolsList = list(symbols.keys())

    for gene in genes.keys():
        species = genes[gene]
        fish = species[3][1] 

        #find lenght of shortest sequence
        minLen = float("inf")
        for i in species:
            if len(i[1]) < minLen:
                minLen = len(i[1])

        #get consensus of mammals for this gene
        profile = [[0]*minLen for i in range(len(symbols))]
        consesus = ""
        for i in range(minLen):
            #build profile one position through all strings
            for k in [0, 1, 2, 4, 5, 6]: #mammals
                s = species[k][1]
                profile[symbols[s[i]]][i] += 1

            #find most common nucleotide
            m = [0, 0]
            for j in range(len(symbols)):
                if profile[j][i] > m[0]:
                    m = [profile[j][i], j]
            consesus += symbolsList[m[1]]

        #distance between consensus and fish is rank
        (optAls1, optAls2) = globalAlingment(fish, consesus)
        eDistance = jukesCantor(optAls1, optAls2)
        ranking.append((gene, eDistance))
    
    #sort decreasing
    ranking = sorted(ranking, key=lambda x: x[1], reverse=True)
    return ranking


if __name__ == "__main__":
    #read file
    f = open("genes.txt", "r")
    data = f.read()
    f.close

    genes = parseGenes(data)

    #rank genes
    start = time.time()
    print("Genes ranking:")
    ranking = rankGenes2(genes)
    [print("\t"+i[0]+" "+str(i[1])) for i in ranking]
    print("Time for ranking: "+str(time.time() - start)+"s.")
    
    #phylogenetic tree for best and worst descriminatory genes
    phylogeneticTree(genes, ranking[0][0]) #best
    phylogeneticTree(genes, ranking[-1][0]) #worst
