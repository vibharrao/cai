import sys
from collections import defaultdict
from statistics import mean
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

'''This program provides functions to calculate codon usage metrics including Relative Synonymous Codon Usage (RSCU), 
the Wi values, and the Codon Adaptation Index (CAI).
Example command line Input: python3 findCAI.py Chiba.fasta NIID123.fasta PRVABC59.fasta >zika.txt'''

class FastAreader:
    def __init__(self, fname=''):
        self.fname = fname
        

    def doOpen(self):
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            header = ''
            sequence = ''
            
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
            header = line[1:].rstrip()
            
            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence

def getCodonTable():
    ''' Retrieve the standard genetic codon table.'''
    codonTable = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    return codonTable

def ParseFasta(filePath):
    '''Parse sequences from a FASTA file.'''
    fastaReader = FastAreader(filePath)
    sequences = {}
    for header, sequence in fastaReader.readFasta():
        sequences[header] = sequence
    return sequences

def calculateRSCU(sequence):
    '''Calculate the Relative Synonymous Codon Usage (RSCU) values for a given sequence.'''
    codonCount = defaultdict(int)  # Initialize a defaultdict to store codon counts
    totalCodons = len(sequence) // 3  # Calculate the total number of codons in the sequence
    codonTable = getCodonTable()
    # Iterate through the sequence in steps of 3 to get each codon
    for i in range(0, totalCodons, 3):
        codon = sequence[i:i+3]  # Extract codon
        if codonTable.get(codon) != '*':  # Exclude stop codons
            codonCount[codon] += 1  # Increment the count for the codon
       
    RSCU = {} # Initialize a dictionary to store RSCU values
    codonCountCopy = dict(codonCount) # Create a copy of codonCount dictionary
    
    # Iterate through each codon and calculate RSCU value
    for codon, count in codonCountCopy.items():
        synonymousCodons = getSynonymousCodons(codon) # Get synonymous codons for the current codon
        if synonymousCodons:  # If there are synonymous codons
            mean_codon_count = mean([codonCount[syn_codon] for syn_codon in synonymousCodons]) # Calculate the mean count of synonymous codons
            RSCU[codon] = count / mean_codon_count  # Calculate RSCU value and store it in the dictionary
        else:  # If there are no synonymous codons
            RSCU[codon] = 0  # Set RSCU value to 0
    return RSCU


def getSynonymousCodons(codon):
    '''Get the synonymous codons for a given codon.'''
    aminoAcid = getCodonTable().get(codon) # Retrieve the amino acid corresponding to the given codon from the codon table
    if aminoAcid is not None:  # Check if the codon corresponds to a valid amino acid
        synonymousCodons = [key for key, value in getCodonTable().items() if value == aminoAcid] # create a list of synonymous codons by iterating over the codon table
        return synonymousCodons  # Return the list of synonymous codons
    else:
        return []  # Return an empty list if the codon is invalid or does not correspond to an amino acid

def calculateWi(RSCUvalues):
    '''Calculate the Wi values for a given set of RSCU values.''' 
    WiValues = {} # Initialize a dictionary to store the Wi values
    RSCUvaluesCopy = dict(RSCUvalues) # Create a copy of the RSCUvalues dictionary
    # Iterate through each codon and its corresponding RSCU value
    for codon, count in RSCUvaluesCopy.items():
        synonymousCodons = getSynonymousCodons(codon)  # Get synonymous codons for the current codon
        meanCodonCount = max([RSCUvalues[synCodon] for synCodon in synonymousCodons if synCodon in RSCUvalues]) # Calculate the maximum RSCU value among synonymous codons
        WiValues[codon] = count / meanCodonCount  # Calculate the Wi value for the current codon and store it in the WiValues dictionary
    return WiValues # Return the dictionary containing the Wi values for each codon


def calculateCAI(RSCUvalues):
    '''Calculate the Codon Adaptation Index (CAI) for a given set of RSCU values.'''
    WiValues = calculateWi(RSCUvalues) # Calculate the Wi values using the provided function
    geometricMean = 1 # Initialize a variable to store the geometric mean
    numCodons = len(WiValues)  # Get the total number of codons
    
    for value in WiValues.values():
        geometricMean *= value # Calculate the geometric mean by multiplying all Wi values

    CAI = geometricMean ** (1 / numCodons) # Calculate CAI by taking the geometric mean to the power of 1/numCodons
    return CAI  # Return the calculated CAI


def calculateCAIref():
    ''''Calculate the reference Codon Adaptation Index (CAI) using predefined RSCU values.'''
    RSCUrefValues = {
    'TTT': 16.9, 'TTC': 20.4, 'TTA': 7.2, 'TTG': 12.6,
    'TCT': 14.6, 'TCC': 17.4, 'TCA': 11.7, 'TCG': 4.5,
    'TAT': 12.0, 'TAC': 15.6, 'TAA': 0.7, 'TAG': 0.5,
    'TGT': 9.9, 'TGC': 12.2, 'TGA': 1.3, 'TGG': 12.8,
    'CTT': 12.8, 'CTC': 19.4, 'CTA': 6.9, 'CTG': 40.3,
    'CCT': 17.3, 'CCC': 20.0, 'CCA': 16.7, 'CCG': 7.0,
    'CAT': 10.4, 'CAC': 14.9, 'CAA': 11.8, 'CAG': 34.6,
    'CGT': 4.7, 'CGC': 10.9, 'CGA': 6.3, 'CGG': 11.9,
    'ATT': 15.7, 'ATC': 21.4, 'ATA': 7.1, 'ATG': 22.3,
    'ACT': 12.8, 'ACC': 19.2, 'ACA': 14.8, 'ACG': 6.2,
    'AAT': 16.7, 'AAC': 19.5, 'AAA': 24.0, 'AAG': 32.9,
    'AGT': 11.9, 'AGC': 19.4, 'AGA': 11.5, 'AGG': 11.4,
    'GTT': 10.9, 'GTC': 14.6, 'GTA': 7.0, 'GTG': 28.9,
    'GCT': 18.6, 'GCC': 28.5, 'GCA': 16.0, 'GCG': 7.6,
    'GAT': 22.3, 'GAC': 26.0, 'GAA': 29.0, 'GAG': 40.8,
    'GGT': 10.8, 'GGC': 22.8, 'GGA': 16.3, 'GGG': 16.4
}

    WirefValues = calculateWi(RSCUrefValues)

    geometricMeanref = 1
    numCodonsref = len(WirefValues)

    for value in WirefValues.values():
        geometricMeanref *= value

    CAIref = geometricMeanref ** (1 / numCodonsref)
    return CAIref

def plotRSCUGrouped(RSCUValues, sampleName, geneName, filename):
    '''Plot RSCU values grouped by amino acids.'''
    codonTable = getCodonTable()
    codons = list(RSCUValues.keys())
    values = list(RSCUValues.values())
    
    # Group codons by amino acid
    aminoAcids = defaultdict(list)
    for codon in codons:
        aminoAcid = codonTable.get(codon)
        if aminoAcid:
            aminoAcids[aminoAcid].append((codon, RSCUValues[codon]))

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))  # Adjust figure size
    bar_width = 0.5  # Adjust the width of each bar
    num_amino_acids = len(aminoAcids)
  
    # Calculate x positions for each amino acid group
    x_positions = np.arange(num_amino_acids)
    
    for i, (aminoAcid, codonValues) in enumerate(aminoAcids.items()):
        for j, (codon, rscu_value) in enumerate(codonValues):
            x_position = x_positions[i] + j * bar_width / len(codonValues)
            ax.bar(x_position, rscu_value, bar_width / len(codonValues), label=codon)
            ax.text(x_position, rscu_value + 0.01, codon, ha='center', va='bottom', fontsize=8)  # Adjust font size
    
    ax.set_xlabel('Amino Acids', fontsize=12)  # Adjust font size
    ax.set_ylabel('RSCU Values', fontsize=12)  # Adjust font size
    ax.set_title(f'RSCU values grouped by amino acids for {geneName} in {sampleName}', fontsize=14)  # Adjust font size
    
    # Set x-axis labels and adjust rotation
    ax.set_xticks(x_positions)
    ax.set_xticklabels(aminoAcids.keys(), rotation=45, fontsize=10)  # Adjust font size and rotation
    
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1), fontsize=8)  # Adjust legend font size
    
    # Adjust margins
    plt.subplots_adjust(bottom=0.15, top=0.85, left=0.1, right=0.9)
    
    plt.savefig(filename)
    plt.close()  # Close the plot to release memory

def compareRSCUAcrossSamples(geneSample1RSCU, geneSample2RSCU, geneSample3RSCU):
    '''Compare RSCU values across multiple gene samples and create a heatmap. '''
    # Create a list of all codons present in the RSCU dictionaries of all samples
    allCodons = list(set(geneSample1RSCU.keys()) | set(geneSample2RSCU.keys()) | set(geneSample3RSCU.keys()))
    codonTable = getCodonTable()

    # Group codons by amino acid
    aminoAcidCodons = defaultdict(list)
    for codon in allCodons:
        aminoAcid = codonTable.get(codon)
        if aminoAcid:
            aminoAcidCodons[aminoAcid].append(codon)

    # Create labels with amino acid names and grouped codons
    labels = []
    for aminoAcid, codons in aminoAcidCodons.items():
        for codon in codons:
            labels.append(f'{codon} ({aminoAcid})')

    # Initialize a matrix to store RSCU values for each codon in each sample
    RSCUMatrix = np.zeros((len(labels), 3))

    for i, codon in enumerate(allCodons):
        RSCUMatrix[i, 0] = geneSample1RSCU.get(codon, 0)
        RSCUMatrix[i, 1] = geneSample2RSCU.get(codon, 0)
        RSCUMatrix[i, 2] = geneSample3RSCU.get(codon, 0)

    # Create a heatmap
    plt.figure(figsize=(12, 8))

    sns.heatmap(RSCUMatrix, cmap="YlGnBu", annot=False, xticklabels=['ZIKV_CHIBA', 'ZIKV_NIID123', 'ZIKV_PRVABC59'], yticklabels=labels)
    plt.xlabel('Gene Samples')
    plt.ylabel('Codons and Amino Acids')
    plt.title('Comparison of RSCU Values Across Gene Samples')
    plt.savefig('heatmap_RSCU_comparison.png')

def main(geneSample1File, geneSample2File, geneSample3File):
    ''' Main function to process gene samples and perform analysis.'''
    geneSample1 = ParseFasta(geneSample1File)
    geneSample2 = ParseFasta(geneSample2File)
    geneSample3 = ParseFasta(geneSample3File)
    
    caiValues = {'ZIKV_CHIBA': [], 'ZIKV_NIID123': [], 'ZIKV_PRVABC59': [], 'Reference': []}
    
    for geneName, geneSequence in geneSample1.items():
        RSCU1 = calculateRSCU(geneSequence)
        CAI1 = calculateCAI(RSCU1)
        plotRSCUGrouped(RSCU1, "Gene Sample 1", geneName, f"{geneName}_RSCU_Gene_Sample_1.png")
        print(f"CAI for {geneName}: {CAI1:.4f}")
        caiValues['ZIKV_CHIBA'].append(CAI1)
    
    for geneName, geneSequence in geneSample2.items():
        RSCU2 = calculateRSCU(geneSequence)
        CAI2 = calculateCAI(RSCU2)
        plotRSCUGrouped(RSCU2, "Gene Sample 2", geneName, f"{geneName}_RSCU_Gene_Sample_2.png")
        print(f"CAI for {geneName}: {CAI2:.4f}")
        caiValues['ZIKV_NIID123'].append(CAI2)
    
    for geneName, geneSequence in geneSample3.items():
        RSCU3 = calculateRSCU(geneSequence)
        CAI3 = calculateCAI(RSCU3)
        plotRSCUGrouped(RSCU3, "Gene Sample 3", geneName, f"{geneName}_RSCU_Gene_Sample_3.png")
        print(f"CAI for {geneName}: {CAI3:.4f}")
        caiValues['ZIKV_PRVABC59'].append(CAI3)
    
    CAIref = calculateCAIref()
    print(f"Reference CAI: {CAIref:.4f}")
    caiValues['Reference'].append(CAIref)

    compareRSCUAcrossSamples(RSCU1, RSCU2, RSCU3)
    
    # Plot CAI values for each gene sample
    fig, ax = plt.subplots()
    sampleNames = list(caiValues.keys())
    caiMeans = [np.mean(values) for values in caiValues.values()]
    ax.bar(sampleNames, caiMeans, width=0.2)  # Adjust the width of the bars as needed
    ax.set_ylabel('Average CAI')
    ax.set_title('Comparison of Average CAI across Gene Samples')
    
    # Save the plot as a PNG file
    plt.savefig('cai_comparison.png')
    
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python program_name.py gene_sample1_file gene_sample2_file gene_sample3_file")
        sys.exit(1)
    
    gene_sample1_file = sys.argv[1]
    gene_sample2_file = sys.argv[2]
    gene_sample3_file = sys.argv[3]
    
    main(gene_sample1_file, gene_sample2_file, gene_sample3_file)
