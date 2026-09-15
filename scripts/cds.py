from Bio import SeqIO

#Read the nucleotide sequence stored in homo_mito fasta file
store = SeqIO.read("data/homo_mito.fasta","fasta")      #store is an object
sequence = store.seq          #sequence is a variable

#CDS coordinates obtained from 1 b) 8527...9207
start = 8527
end = 9207

#Extract CDS from the genome
cds = sequence[start - 1:end]  #The end point was not subtracted by 1 because 9207 will not be included onlt till 9206 
#and in python indexing is from 0 so we need till 9206

#Result
print("The CDS sequence is:")
print(cds)

print("The length of CDS is:", len(cds),"nucleotides")
print("The remainder when length is divided by 3 is:", len(cds) % 3)

if len(cds)%3 ==0:
   print("Length of CDS is divisible by 3.")
else:
   print("The CDS length is not divisible by 3.")

#Dictionary for vertebrate mitochondrial genetic code as listed by NCBI
mito_genetic_code = {  "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "", "TAG": "",
    "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",

    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",

    "ATT": "I", "ATC": "I", "ATA": "M", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "", "AGG": "",

    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

# "" represents a stop codon
#Translating the sequence
protein = ""
for i in range(0, len(cds), 3):
    codon = str(cds[i:i+3])
    amino_acid = mito_genetic_code[codon]
    protein += amino_acid

#Results
print("First 30 translated residues are:", protein[:30])
print("Total length of the protein is",len(protein))

#Save translated protein to a new fasta file
with open("data/mito_residue.fasta","w") as file:
       file.write(">MT ATP Protein\n")
       file.write(str(protein) + "\n")

#Read the translated protein sequence
with open("data/mito_residue.fasta") as sh:
#Combining all protein sequence lines into one continuous sequence without space
  my_trans_seq = "".join(sh.read().splitlines()[1:])
print("The sequence translated now is:\n" + my_trans_seq + "\n")

#Read the protein file downloaded in 1 c)
with open("data/Protein_ATP_synthase_F0_a.fasta") as sh:
  deposited_seq = "".join(sh.read().splitlines()[1:])
print("The already available sequence is:\n" + deposited_seq)

#Comparing both
difference = False
for i in range(len(my_trans_seq)):
  if my_trans_seq[i] != deposited_seq[i]:
        print("Difference at position:", i+1)
        difference = True

if difference == False:
        print("Both sequences are identical.")

#Cross-check using Biopython
biopython_protein = str(cds.translate(table=2))  #cds is the variable which has extracted CDS region
#table 2 is taken because here human mitochondrial gene is being studied
#str used to connvert the type that Biopython returns into
print("Biopython translation:\n" + biopython_protein + "\n")

print("My translated protein:\n" + my_trans_seq + "\n")

#Compare my translated sequence with Biopython translation
if my_trans_seq == biopython_protein:
    print("Both translations match.")
else:
    print("The translations are different from each other.")

#Length of translated and deposited protein
print("The length of my translated protein is:")
print(len(my_trans_seq))
print("\n" + "The length of the deposited protein is:")
print(len(deposited_seq))
print("\n" + "Difference between both the protein sequences is:")
print(len(my_trans_seq) - len(deposited_seq))

# For pipeline file
with open("results/MT_ATP6_results.txt", "w") as fh:
    fh.write("MT ATP6\n")
    fh.write("Category B - Mitochondrial\n")
    fh.write(str(len(cds)) + "\n")
    fh.write(str(len(my_trans_seq)) + "\n")
    fh.write(str(len(deposited_seq)) + "\n")
    fh.write("Y\n")
    fh.write("N/A\n")

print("MT_ATP6_results.txt created")
def show_summary(files):

    print(f"{'Gene':<12}{'Category':<30}{'CDS Length':<15}"
          f"{'Translated':<15}{'Deposited':<15}"
          f"{'Identical':<12}{'First Mismatch'}")

    for filename in files:

        with open(filename) as f:
            data = [line.strip() for line in f]

        print(f"{data[0]:<12}{data[1]:<30}{data[2]:<15}"
              f"{data[3]:<15}{data[4]:<15}{data[5]:<12}{data[6]}")
show_summary([
    "results/MT_ATP6_results.txt"])
