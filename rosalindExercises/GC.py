#read file
f = open("10data.txt", "r")
data = f.read()
f.close

#split FASTA format
fasta = data.split(">")[1:]
dnaDict = {}

#anaylze
maxLabel = ""
maxCG = 0
for i in fasta:
    #dnaDict[s.split("\n", 1)[0]] = s.split("\n", 1)[1]
    label = i[:14]
    strand = i[14:]
    strand = strand.replace('\n', '')

    n = len(strand)
    cg = 0

    for nd in strand:
        if nd == "C" or nd == "G":
            cg += 1

    cg = (cg/n)*100
    if cg > maxCG:
        maxCG = cg
        maxLabel = label

print("%s%f" %(maxLabel, maxCG))



