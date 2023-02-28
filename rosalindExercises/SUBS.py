f = open("12data.txt", "r")
data = f.read()
f.close()
(s, t) = data.split("\n")

n = len(t)
occurences=""

#s[:-(n-1)]...n-1 zato ker nocemo izbrisati zadnjih n ampak "enega" pustiti da lahko zajamemo celotni zadnji podniz
for i in range(len(s[:-(n-1)])):
    if s[i:i+n] == t:
        #i+1 ker v printu hocemo steti od 1 ne od 0
        occurences += str(i+1)+" "

print(occurences)