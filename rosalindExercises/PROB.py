import math

#read and parse data
f = open("16data.txt", "r")
data = f.read()
(s, GCs) = data.split("\n")
GCs = GCs.split(" ")
GCs = [float(i) for i in GCs]

B = [] #result array

for gc in GCs:
    #get probabilities for constructing random DNA string
    probBp = {"G": gc/2, "C": gc/2, "A":(1-gc)/2, "T":(1-gc)/2}

    #go throug s and compute log probability for such a string occuring randomly
    prob = 0
    for i in s:
        prob += math.log10(probBp[i])

    B.append(round(prob, 3))

#write solution to file
f = open("16result.txt", "w")
[f.write(str(i)+" ") for i in B]
f.close