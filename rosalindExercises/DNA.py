strand = "AGATTAACGCGAGTGATAGTCTCAAGGTTCGCACCGACTGCGCTTGTAGCCCCAACTGGACTGTAGGTGATGGCACATAGTCCCTAAACATCCGTCTAACCCCTTTGGCCTTAGGCTGCACGTTCAGGGTAGTTATACCATCCGGCACAACCCCCCAGCTGACACCGAAAAGCGAGTGGGACGGATTGAAATCCAGATGTCCGACTCGTTGATCATTTACCTGGTTCGAATGAATAACCCTGGACTATACTATCTAGTTTGATCTCGCTCAGGACCCGGTAGATGACATAAGGAGGGCTGAAGTCACACGGCTTCTCACCGTGAATAGACTACTAAAGACACATCAGTAGGCTGGAGCCCTATATGTGGGACCCGGTTGGATCGCGACTTCGTTATTCGAGCTGGCAACTATCGAAAAGTGGTGGGAGTTGGGCCCGCTATCCGCATCCAAGTAATTGATCCTGGCGGAACTCGCCGCTACCAGTGGGCTTCGTCATACGGTTAGAGCTCTCGTTAGAGGTAGAATCTGAGAGGGGCTACGAGATAGATTCAGCAGTCACGACAGCTCTTCCAAAGAAGCGACTGCCAGTAAATTCGGGACAGCCTAATCGCCGGTGCCGGGTGGTGAAGGCTTACCGCATGTGAGTTGTTTGTGCTGGCAAATTCAGTGTCAAAGGGCAACACCGTGAGTAATCGATGGGAGCCCCACCCGACACACTCGACGAAGTGGCAATCCGTGGCGGAGGGATGTGAGCGACTTCGAAAGTTAGAGCCAAAGATAATAAGTGACCAACGGTGGGTCTTGAGAAAGCTCTTCCCGTGCCTATCGCGGATCCGCTGCGGCATCTTTTCAGGAGCAACACCGTAGC"
a, c, g, t = 0,0,0,0

for nt in strand:
    if nt == "A":
        a += 1
    elif nt == "C":
        c += 1
    elif nt == "G":
        g += 1
    else:
        t += 1
    
print("%d %d %d %d" %(a,c,g,t))