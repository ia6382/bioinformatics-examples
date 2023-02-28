def createKmers(alphabet, prefix, k):
    if k == 0:
        return [prefix]
    else:
        kmers = []
        for i in range(len(alphabet)):
            newPrefix = prefix + alphabet[i]
            kmers = kmers + createKmers(alphabet, newPrefix, k-1)
        return kmers

#read data
f = open("14data.txt", "r")
data = f.read()

#split string
alphabet, k = data.split("\n")
alphabet = alphabet.split(" ")

#create all posible kmers
kmers = createKmers(alphabet, "", int(k))

#write solution to file
f = open("14result.txt", "w")
[f.write(i+"\n") for i in kmers]
f.close