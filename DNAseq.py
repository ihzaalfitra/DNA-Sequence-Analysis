from Bio import SeqIO

def ReadFasta(file_path):
    """Reads a FASTA file and returns the DNA sequence as a string."""
    with open(file_path, "r") as fasta_file:
        record = next(SeqIO.parse(fasta_file, "fasta"))  # Read first sequence
        return str(record.seq)  # Convert sequence to string

def PatternCount(Text, Pattern):
    """Finds the frequency of pattern within text"""
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):  # Sliding window
        # print(Text[i:i+len(Pattern)])
        if Text[i:i+len(Pattern)] == Pattern:
            count += 1
    return count

def FrequentWords(Text, k):
    """Finds the most frequent k-mers in Text."""
    FrequentPatterns = set()
    Count = []

    # Count occurrences of each k-mer
    for i in range(len(Text) - k + 1):
        Pattern = Text[i:i+k]
        Count.append(PatternCount(Text, Pattern))

    maxCount = max(Count)

    # Find all k-mers with maxCount
    for i in range(len(Text) - k + 1):
        if Count[i] == maxCount:
            FrequentPatterns.add(Text[i:i+k])  # Add to set (removes duplicates)

    return FrequentPatterns

def FrequencyTable(Text, k):
    """Generates a frequency table (dictionary) of all k-mers in the given Text."""
    freqMap = {}  # Empty dictionary to store k-mer counts
    n = len(Text)  # Length of the text

    # Iterate over all k-mers in the Text
    for i in range(n - k + 1):
        Pattern = Text[i:i + k]  # Extract k-mer

        if Pattern not in freqMap:  # If Pattern not in dictionary, initialize it
            freqMap[Pattern] = 1
        else:  # Otherwise, increment count
            freqMap[Pattern] += 1

    return freqMap  # Return dictionary of k-mer frequencies

def MaxMap(freqMap):
    """Finds the maximum frequency value in the frequency map."""
    return max(freqMap.values())  # Returns the highest k-mer frequency


def BetterFrequentWords(Text, k):
    """Finds the most frequent k-mers in Text using a more efficient approach."""
    FrequentPatterns = []  # List to store the most frequent k-mers
    freqMap = FrequencyTable(Text, k)  # Compute k-mer frequency table
    maxCount = MaxMap(freqMap)  # Find the maximum k-mer count

    # Find all k-mers that have the max count
    for Pattern, count in freqMap.items():
        if count == maxCount:
            FrequentPatterns.append(Pattern)

    return FrequentPatterns  # Return the list of most frequent k-mers
def ComplementaryBase(Base):
    complements = {
        'A': 'T',  # Adenine pairs with Thymine
        'T': 'A',  # Thymine pairs with Adenine
        'C': 'G',  # Cytosine pairs with Guanine
        'G': 'C'   # Guanine pairs with Cytosine
    }
    return complements.get(Base.upper())

def RevereComplement(Pattern):
    ReverseComplement=""
    for i in range(len(Pattern)):
        ReverseComplement = ComplementaryBase(Pattern[i]) + ReverseComplement
    return ReverseComplement

def PatternPositions(Pattern, Genome):
    """Finds the frequency of pattern within text"""
    positions = []

    for i in range(len(Genome) - len(Pattern) + 1):  # Sliding window
        # print(Text[i:i+len(Pattern)])
        if Genome[i:i+len(Pattern)] == Pattern:
            positions.append(i)
    return positions

# Example Usage
# Text = "ATGATGATGCTG"
# k = 3
# freqMap = FrequencyTable(Text, k)

# Print the frequency table
# for pattern, count in freqMap.items():
#     print(f"{pattern}: {count}")

# Example Usage
Vibrio_cholerae = ReadFasta("Vibrio_cholerae.fasta")
Genome = Vibrio_cholerae
Pattern = "ATGATCAAG"
k = 14

# print(BetterFrequentWords(Text, k))
# print(Text)
print(PatternPositions(Pattern, Genome))
# print(RevereComplement(Pattern))
# print(Text[1:3])
# print(range(len(Text) - len(Pattern) + 1))

"""
Notes:
Synthesis of a complementary strand on a template strand
See Meselson and Stahl in 1958

The beginning and end of a DNA strand are denoted 5’ (pronounced “five prime”) and 3’ (pronounced “three prime”)

Ori-finder
https://tubic.org/doric/tools


https://www.nature.com/articles/171737a0

book companion:
https://www.cl.cam.ac.uk/~pl219/Bioinformatics2015.pdf
courses:
https://www.coursera.org/specializations/bioinformatics
https://cogniterra.org/lesson/29849/step/1?unit=21946
"""
