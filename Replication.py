def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text) - (len(Pattern)-1)):
        if Text[i:(i+len(Pattern))] == Pattern:
            count += 1
    return count


def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = PatternCount(Text, Pattern)
    return freq


def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words


def ReverseComplement(Pattern):
    rev_seq = Reverse(Pattern)
    comp_seq = Complement(rev_seq)
    return comp_seq


def Reverse(Pattern):
    return Pattern[::-1]


def Complement(Pattern):
    seq = []
    for base in Pattern:
        if base == "A":
            seq.append("T")
        elif base == "T":
            seq.append("A")
        elif base == "C":
            seq.append("G")
        elif base == "G":
            seq.append("C")
    seq = "".join(seq)
    return seq


def PatternMatching(Pattern, Genome):
    positions = []
    k = len(Pattern)
    n = len(Genome)
    for i in range(n-k+1):
        if Pattern == Genome[i:i+k]:
            positions.append(i)
    return positions


def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(ExtendedGenome[i:i+(n//2)], symbol)
    return array


def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    array[0] = PatternCount(Genome[0:n//2], symbol)
    for i in range(1, n):
        array[i] = array[i-1]
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array


def SkewArray(Genome):
    skew = [0]
    for i in range(len(Genome)):
        if Genome[i] == "C":
            skew.append(skew[i]-1)
        elif Genome[i] == "G":
            skew.append(skew[i]+1)
        else:
            skew.append(skew[i])
    return skew


def MinimumSkew(Genome):
    min_positions = []
    skew = SkewArray(Genome)
    skew_min = min(skew)
    for i in range(len(Genome)):
        if skew[i] == skew_min:
            min_positions.append(i)
    return min_positions


def HammingDistance(p, q):
    difference = 0
    if len(p) == len(q):
        for i in range(len(p)):
            if p[i] != q[i]:
                difference += 1
        return difference
    else:
        print("Not equal length")

def ApproximatePatternMatching(Text, Pattern, d):
    positions = []
    k = len(Pattern)
    n = len(Text)
    for i in range(n-k+1):
        if HammingDistance(Pattern, Text[i:i+k]) <= d:
            positions.append(i)
    return positions


