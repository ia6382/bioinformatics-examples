def hDistance(s1, s2):
    hammingDistance = 0
    for (i, j) in zip(s1, s2):
        if i != j:
            hammingDistance += 1
    
    return hammingDistance

if __name__ == "__main__":

    #read file
    f = open("20data.txt", "r")
    data = f.read()
    f.close

    (s1, s2) = data.split("\n")
    hd = hDistance(s1, s2)

    print(hd)
