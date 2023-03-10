strand = "CTCCGACCCATTAATCTAGGCTTGAACTCAATTGTCTCCATGACGCAGCCGTGAGTTGAATCCAGCATGATGAAGTATCGTGCTGGTAGGGTGTCAGCAACAGGTTTCACCAAATACACGGCCTTAGGCGGACTCGTTCGTCGGTGAGGATACTACTTGACGCGGGCTGCTTTACCTGCCATCTGAGACCAGACTATACACTGAGCGGCTTGTGCTTTCGCTGGCGAGTGGATATGTGGGCCCGTCAACTGTGTTCGAGCCCAGCGTCTGTGGCTACCACCGATGCCACACTCTAAGTCGTGCACATTCGGTGCGGAAGCGCGTGCTGCGTCGCCGAGCTACAGACTATTCCTCCATAAGGCTCAACAGTCTTTTATAGAGTTCGAGTTAAACCAGAGACGAGGCAACTAGACGGCAGCTACAGAGGGTACAGATATGGTTGATAGCGCTGGTCCCACGAGATACCACGTCTGACCAAAAGAGAACTGAGCCGCGAGTCTTTGAGATTTAATACCATTTCTATTGTTTCTGAAACACACTGTTTTAAAACAGCACCGAATATTCCCAACAGTTACAATAGTGTGGGGTCCGTATGCACTCATAGTTCTCCATATCCAGGGACGGAGTCATCTTCCGCGCCTATTGGAGTAGCCACACGTTTTACCGAGCGTTGATACATACAGTGGTTCACTGCTGGTTATGAGTGCATCCGAAAGCTCGCCTACCGTTGTGATTGGCTATGTTCATACCACTACTGTGACTTTGGCCCTGTGAGTACCAATTGGAGCTTTCGACCTCGGATAGTAATTTCACCACATGGTCATGGCTAGTGCGGCGCTTTCATTCAACATTCCTCGC"
cStrand = ""

for nt in strand:
    if nt == "A":
        cStrand += "T"
    elif nt == "C":
        cStrand += "G"
    elif nt == "G":
        cStrand += "C"
    else:
        cStrand += "A"

cStrand = cStrand[::-1]
print(cStrand)