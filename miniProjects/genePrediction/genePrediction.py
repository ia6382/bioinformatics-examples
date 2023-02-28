from Bio import Entrez
from Bio import SeqIO
import matplotlib.pyplot as plt
from statistics import median

dnaCodonTable ={
"TTT":"F","CTT":"L","ATT":"I","GTT":"V",
"TTC":"F","CTC":"L","ATC":"I","GTC":"V",
"TTA":"L","CTA":"L","ATA":"I","GTA":"V",
"TTG":"L","CTG":"L","ATG":"M","GTG":"V",
"TCT":"S","CCT":"P","ACT":"T","GCT":"A",
"TCC":"S","CCC":"P","ACC":"T","GCC":"A",
"TCA":"S","CCA":"P","ACA":"T","GCA":"A",
"TCG":"S","CCG":"P","ACG":"T","GCG":"A",
"TAT":"Y","CAT":"H","AAT":"N","GAT":"D",
"TAC":"Y","CAC":"H","AAC":"N","GAC":"D",
"TAA":"s*","CAA":"Q","AAA":"K","GAA":"E",
"TAG":"s*","CAG":"Q","AAG":"K","GAG":"E",
"TGT":"C","CGT":"R","AGT":"S","GGT":"G",
"TGC":"C","CGC":"R","AGC":"S","GGC":"G",
"TGA":"W","CGA":"R","AGA":"R","GGA":"G",
"TGG":"W","CGG":"R","AGG":"R","GGG":"G"
}

def downloadData():
    Entrez.email = "ia6382@student.uni-lj.si"
    id = "NC_000908"
    with Entrez.efetch(db="nucleotide", id=id, rettype="fasta", retmode="text") as handle:
        data = SeqIO.read(handle, "fasta")

    return data

def downloadAnnotations():
    id = "NC_000908"
    with Entrez.efetch(db="nucleotide", rettype="gb", id=id) as handle:
        rec = SeqIO.read(handle, "gb")

    return rec

def translate(dna):
    acids = []
    for i in range(0, len(dna), 3):
        codon = dna[i:i+3]
        aminoAcid = dnaCodonTable[codon]
        acids.append(aminoAcid)

    acids = ''.join(acids)
    return acids

def reverseComplement(dna):
    dnaC = []
    for nt in dna[::-1]:
        if nt == "A":
            dnaC.append("T")
        elif nt == "C":
            dnaC.append("G")
        elif nt == "G":
            dnaC.append("C")
        else:
            dnaC.append("A")

    dnaC = ''.join(dnaC)    
    return dnaC

def findProteins(dna, dnaRC):
    #go through all possitions: 3 for different start nt *2 for complement
    proteins = set() #all possible proteins
    for s in (dna, dnaRC):
        for frameStart in range(3):
            length = len(s) - (len(s)-frameStart)%3
            cutS = s[frameStart:length]

            #split tranlated s at stop codons (hack: stop codod marked as *s)
            acids = translate(cutS)
            potentialProteins = acids.split("*")
            
            #cut of acids before start acid M (=codon AUG)
            for p in potentialProteins:
                if p[-1] != "s": #to make sure it has a stop codon at the end
                    continue
                for (j, a) in enumerate(p):
                    if a == "M":
                        proteins.add(p[j:-1]) #ignore the last letter "s"
                        break #only save the longest one => much better precision, but not all possible genes are found
    return proteins

def statisticsOfAllGenes(genom, sequences):
    lenghts = []
    for i in sequences:
        if i.type == "gene":
            gene = str(i.extract(genom).seq.translate()).replace("*", "")
            lenghts.append(len(gene))

    print("\nall genes in data: "+ str(len(lenghts)))
    print("lenght of longest gene (in codons): "+str(max(lenghts)))
    print("lenght of shortest gene (in codons): "+str(min(lenghts)))
    print("median lenght of genes (in codons): "+str(median(lenghts)))


if __name__ == "__main__":
    print("downloading genome ...")
    data = downloadData()

    print("downloading annotations ...")
    rec = downloadAnnotations()

    print("analyzing ...")
    dna = data.seq._data
    dnaRC = reverseComplement(dna)
    
    statisticsOfAllGenes(data, rec.features)

    #create set of all annotated proteins
    realProteins = set()
    lenghts = []
    for i in rec.features:
        if i.type == 'CDS':
            if "translation" in i.qualifiers:
                protein = i.qualifiers["translation"][0]
                realProteins.add(protein)
                lenghts.append(len(protein))

    #report statistics
    print("\ngenes in data that translate to proteins: "+str(len(realProteins)))
    print("lenght of longest gene (in codons): "+str(max(lenghts)))
    print("lenght of shortest gene (in codons): "+str(min(lenghts)))
    print("median lenght of genes (in codons): "+str(median(lenghts)))

    #find all possible proteins
    myProteins = findProteins(dna, dnaRC)
    
    #measure the model effectivness in relation to min length L of protein
    recall = []
    precision = []
    numberOfProteins = []
    r = range(37, 1805)
    for L in r:
        myProteinsL = set(filter(lambda x: len(x) >= L, myProteins))
        
        fp = len(myProteinsL - realProteins)
        tp = len(myProteinsL & realProteins)
        fn = len(realProteins - myProteinsL)

        recall.append(tp/(tp+fn))
        precision.append(tp/(tp+fp))
        numberOfProteins.append(len(myProteinsL))

        if L == 50:
            print("\nL = %d. Recall: %.3f, Precision: %.3f"  %(L, tp/(tp+fn), tp/(tp+fp)))
        elif L == 125:
            print("L = %d. Recall: %.3f, Precision: %.3f"  %(L, tp/(tp+fn), tp/(tp+fp)))

    #plot graph
    x = list(r)
    plt.plot(x, recall, "-b", label='občutljivost')
    plt.plot(x, precision, "-r", label='natančnost')
    plt.xlabel("L")
    plt.title("Uspešnost iskanja genov v genomu bakterije Mycoplasma genitalium ")
    plt.legend(loc="upper left")
    plt.show()

    plt.plot(x, numberOfProteins, "-g", label='number of found genes')
    plt.xlabel("L")
    plt.title("Number of found genes")
    plt.ylabel("#genes")
    plt.show()