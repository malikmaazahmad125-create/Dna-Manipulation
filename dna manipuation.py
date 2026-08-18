import numpy as np
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns


print("." * 30)
print("TITLE: DNA MANIPULATION")
print("." * 30)


# ==========================================
# DNA SEQUENCE
# ==========================================

dna_sequence = "ATGCGATCG"

length = len(dna_sequence)

print("length of dna_seq. :", length)

seq = dna_sequence[0:4]

print("first four characters :", seq)


# GC CONTENT

count = dna_sequence.count("G") + dna_sequence.count("C")

gc_content = (count / length) * 100

print("GC_CONTENT:", gc_content, "%")


# ==========================================
# DNA COMPLEMENT
# ==========================================

def DNA_complement(seq):

    complement = {

        "A": "T",

        "T": "A",

        "G": "C",

        "C": "G"

    }

    reverse_complement = seq[::-1]

    rev_com = ''

    for base in reverse_complement:

        rev_com = rev_com + complement[base]

    return rev_com


result = DNA_complement(dna_sequence)

print("reverse of DNA_COMPLEMENT IS :", result)


# ==========================================
# CODON ANALYSIS
# ==========================================

print("." * 15, "codon analysis", "." * 15)


def count_codon(seq):

    count = 0

    for i in range(0, len(seq) - 2, 3):

        count += 1

    return count


count_result = count_codon(dna_sequence)

print("number of codons is :", count_result)


# ==========================================
# FILE CREATION
# ==========================================

print("." * 15, "file creation", "." * 15)


with open("sample.fasta", "w") as f:

    f.write(
        ">gene1\n"
        "TGCTAGCT\n"
        ">gene2\n"
        "CCTTGGAA\n"
    )


with open("sample.fasta", "r") as f:

    fasta_content = f.read()


print("FASTA FILE CONTENT IS :\n", fasta_content)


# ==========================================
# NUMPY STATISTICS FUNCTIONS
# ==========================================

print("." * 15, "NUMPY STATISTICS FUNCTIONS", "." * 15)


dna_array = np.array(list(dna_sequence))

print("dna as numpy array is:", dna_array)


encoding = {

    "A": 0,

    "T": 1,

    "G": 2,

    "C": 3

}


encoding_seq = np.array(

    [encoding[base] for base in dna_sequence]

)

print("encoding seq. is :", encoding_seq)


gene_expression = np.array(

    [12.5, 15.3, 9.8, 22.1, 18.4]

)


print("the mean of gene_expression is :", np.mean(gene_expression))

print("the median of gene_expression is :", np.median(gene_expression))

print("the standard deviation of gene_expression is:", np.std(gene_expression))

print("the maximum number gene_expression is :", np.max(gene_expression))


# ==========================================
# PANDAS DATAFRAME AND DATA CLEANING
# ==========================================

print("." * 15, "PANDAS DATAFRAME AND DATA CLEANING", "." * 13)


data = pd.DataFrame({

    "Gene": ["BRCA1", "TP53", "EGFR", "MYC"],

    "Expression": [12.5, None, 22.1, 18.4],

    "Condition": ["Control", "Treatment", "Control", "Treatment"]

})


print("dataframe is :")

print(data)


data["Expression"] = data["Expression"].fillna(

    data["Expression"].mean()

)


print("TOPIC 7: Cleaned DataFrame:\n")

print(data)


# ==========================================
# BENEFITS OF DESCRIBE()
# ==========================================

print("." * 13, "BENEFITS OF COMMAND DESCRIBE() IN PANDAS", "." * 13)


data_expression = pd.Series(

    [45, 70, 80, 50, 90]

)

print(data_expression.describe())


print(list(range(3)))


# ==========================================
# DATA VISUALIZATION
# MATPLOTLIB + SEABORN
# ==========================================

print("." * 15, "DATA VISUALIZATION", "." * 15)


sns.set_theme(style="whitegrid")


# ==========================================
# 1. DNA NUCLEOTIDE COUNT
# SEABORN BAR CHART
# ==========================================

nucleotide_counts = Counter(dna_sequence)

nucleotides = ["A", "T", "G", "C"]

counts = [

    nucleotide_counts[n]

    for n in nucleotides

]


plt.figure(figsize=(8, 5))

sns.barplot(

    x=nucleotides,

    y=counts

)

plt.title("DNA Nucleotide Composition")

plt.xlabel("Nucleotides")

plt.ylabel("Count")

plt.show()


# ==========================================
# 2. DNA NUCLEOTIDE PERCENTAGE
# MATPLOTLIB PIE CHART
# ==========================================

percentages = [

    (nucleotide_counts[n] / length) * 100

    for n in nucleotides

]


plt.figure(figsize=(7, 7))

plt.pie(

    percentages,

    labels=nucleotides,

    autopct="%1.1f%%",

    startangle=90

)

plt.title("DNA Nucleotide Percentage")

plt.show()


# ==========================================
# 3. GENE EXPRESSION
# SEABORN BAR CHART
# ==========================================

samples = [

    "Sample 1",

    "Sample 2",

    "Sample 3",

    "Sample 4",

    "Sample 5"

]


plt.figure(figsize=(9, 5))

sns.barplot(

    x=samples,

    y=gene_expression

)

plt.title("Gene Expression Levels")

plt.xlabel("Samples")

plt.ylabel("Expression Level")

plt.show()


# ==========================================
# 4. GENE EXPRESSION TREND
# MATPLOTLIB LINE GRAPH
# ==========================================

plt.figure(figsize=(9, 5))

plt.plot(

    samples,

    gene_expression,

    marker="o"

)

plt.title("Gene Expression Trend")

plt.xlabel("Samples")

plt.ylabel("Expression Level")

plt.grid()

plt.show()


# ==========================================
# 5. DOSE? NO — GENE DATA DISTRIBUTION
# SEABORN HISTOGRAM
# ==========================================

plt.figure(figsize=(8, 5))

sns.histplot(

    gene_expression,

    bins=5,

    kde=True

)

plt.title("Gene Expression Distribution")

plt.xlabel("Expression Level")

plt.ylabel("Frequency")

plt.show()


# ==========================================
# 6. PANDAS DATA VISUALIZATION
# GENE EXPRESSION BY GENE
# ==========================================

plt.figure(figsize=(8, 5))

sns.barplot(

    data=data,

    x="Gene",

    y="Expression",

    hue="Condition"

)

plt.title("Gene Expression by Gene")

plt.xlabel("Gene")

plt.ylabel("Expression")

plt.show()


print("\n" + "=" * 40)

print("DNA MANIPULATION PROJECT COMPLETED!")

print("=" * 40)
