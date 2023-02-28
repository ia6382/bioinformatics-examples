#read file
f = open("13data.txt", "r")
data = f.read()

#split FASTA format
fasta = data.split(">")[1:]
fasta = fasta[0]
label = fasta[:14]
strand = fasta[14:]
strand = strand.replace('\n', '')

#build 4-mers
kmers = {}
alphabet = "ATCG"
for i in alphabet:
    for j in alphabet:
        for k in alphabet:
            for l in alphabet:
                kmers[i+j+k+l] = 0
k = 4

#analyze strand
for i in range(len(strand[:-(k-1)])):
    kmers[strand[i:i+k]] += 1

#write lexicographically
f = open("13result.txt", "w")

for i,(k, v) in enumerate(sorted(kmers.items(), key=lambda x: x[0])):
    f.write(str(v)+" ")

f.close