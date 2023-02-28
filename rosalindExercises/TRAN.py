def parseFasta(data):
    #split FASTA format
    fasta = data.split(">")[1:]
    S = []
    for i in fasta:
        label = i[:14]
        strand = i[14:]
        strand = strand.replace('\n', '')
        S.append(strand)

    return S

if __name__ == "__main__":
    #read file
    f = open("39data.txt", "r")
    data = f.read()
    f.close

    (s1, s2) = parseFasta(data)

    transitions = 0
    transversions = 0
    for n1, n2 in zip(s1, s2):
        if n1 != n2:
            if (n1 == "A" and n2 == "G") or (n1 == "G" and n2 == "A") or (n1 == "C" and n2 == "T") or (n1 == "T" and n2 == "C"):
                transitions += 1
            elif (n1 == "A" and n2 == "C") or (n1 == "C" and n2 == "A") or (n1 == "G" and n2 == "T") or (n1 == "T" and n2 == "G") or (n1 == "A" and n2 == "T") or (n1 == "T" and n2 == "A") or (n1 == "G" and n2 == "C") or (n1 == "C" and n2 == "G"):
                transversions += 1

    Rtt = transitions/transversions
    print(Rtt)