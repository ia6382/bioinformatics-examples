f = open("15data.txt", "r")
data = f.read()

#split FASTA format
fasta = data.split(">")[1:]
fasta = fasta[0]
label = fasta[:14]
strand = fasta[14:]
strand = strand.replace('\n', '')

#KMP algoritem, samo Fail tabela, ki vsebuje dolzine prefixov
T = [0]*len(strand)
T[0] = 0 # taksan je navada, v wiki je -1 in je alg drugacen

pos = 1 #indeks v T
cnd = 0 #indeks naslednjega kandidata zaporedja v string

while pos < len(T):
    if strand[pos] == strand[cnd]:
        cnd = cnd + 1  #prefix se nadaljuje povecaj cnd 
    else:
        while cnd > 0 and strand[pos] != strand[cnd]: #koliko nazaj moremo, da so spet enaki
            cnd = T[cnd-1] # vrednost je indeks konca do sedaj najdaljsega zaporedja

        if strand[pos] == strand[cnd]: # ce sta enaka zacni nov prefix
            cnd = cnd + 1
    
    #zapisi dolzino prefixa
    T[pos] = cnd
    #nova iteracija
    pos = pos + 1


#write solution to file
f = open("15result.txt", "w")
[f.write(str(i)+" ") for i in T]
f.close