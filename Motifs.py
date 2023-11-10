def Count(*Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count


def Profile(*Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    matrix = Count(*Motifs)
    for symbol in "ACGT":
        profile[symbol] = []
        for i in range(k):
            profile[symbol].append(matrix[symbol][i]/t)
    return profile


def Consensus(*Motifs):
    k = len(Motifs[0])
    count = Count(*Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus


def Score(*Motifs):
    score = 0
    t = len(Motifs)
    k = len(Motifs[0])
    count = Count(*Motifs)
    consensus = Consensus(*Motifs)
    for j in range(k):
        score += (t - count[consensus[j]][j])
    return score


def Pr(Text, Profile):
    p = 1
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]
    return p


def MaxProbKmers(text, k, profile):
    n = len(text)
    probability = []
    for i in range(n-k+1):
        probability.append(Pr(text[i:i+k], profile))
    return max(probability)


def ProfileMostProbableKmer(text, k, profile):
    n = len(text)
    probability = []
    for i in range(n-k+1):
        probability.append(Pr(text[i:i+k], profile))
    m = max(probability)
    for i in range (n-k+1):
        if Pr(text[i:i+k], profile) == m:
            return text[i:i+k]


def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(*Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(*Motifs) < Score(*BestMotifs):
            BestMotifs = Motifs
    return BestMotifs


#Laplace’s Rule of Succession (avoid zero count)
def CountWithPseudocounts(*Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = {}
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(1)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count


def ProfileWithPseudocounts(*Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    matrix = CountWithPseudocounts(*Motifs)
    for symbol in "ACGT":
        profile[symbol] = []
        for i in range(k):
            profile[symbol].append(matrix[symbol][i]/(t+4))
    return profile


def ConsensusWithPseudocounts(*Motifs):
    k = len(Motifs[0])
    count = CountWithPseudocounts(*Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

def ScoreWithPseudocounts(*Motifs):
    score = 0
    t = len(Motifs)
    k = len(Motifs[0])
    count = CountWithPseudocounts(*Motifs)
    consensus = ConsensusWithPseudocounts(*Motifs)
    for j in range(k):
        score += (t - count[consensus[j]][j])
    return score

def PrWithPseudocounts(Text, Profile):
    p = 1
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]
    return p

def ProfileMostProbableKmerWithPseudocounts(text, k, profile):
    n = len(text)
    probability = []
    for i in range(n-k+1):
        probability.append(PrWithPseudocounts(text[i:i+k], profile))
    m = max(probability)
    for i in range (n-k+1):
        if PrWithPseudocounts(text[i:i+k], profile) == m:
            return text[i:i+k]

def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(*Motifs[0:j])
            Motifs.append(ProfileMostProbableKmerWithPseudocounts(Dna[j], k, P))
        if ScoreWithPseudocounts(*Motifs) < ScoreWithPseudocounts(*BestMotifs):
            BestMotifs = Motifs
    return BestMotifs


def Motifs(Profile, Dna):
    motifs = []
    n = len(Dna)
    k = len(Profile['A'])
    for i in range(n):
        motifs.append(ProfileMostProbableKmer(Dna[i], k, Profile))
    return motifs


import random


def RandomMotifs(Dna, k, t):
    rand_motifs = []
    for i in range(t):
        r = random.randint(0, len(Dna[0])-k)
        rand_motifs.append(Dna[i][r:r+k])
    return rand_motifs


def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(*M)
        M = Motifs(Profile, Dna)
        if Score(*M) < Score(*BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs


def Normalize(Probabilities):
    P_sum = sum(Probabilities.values())
    normalized = {key: value / P_sum for key, value in Probabilities.items()}
    return normalized


def WeightedDie(Probabilities):
    p = random.uniform(0, 1)
    cumulative_prob = 0
    for kmer, prob in Probabilities.items():
        cumulative_prob += prob
        if p < cumulative_prob:
            return kmer


def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {}
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)


# Input:  Integers k, t, and N, followed by a collection of strings Dna
# Output: GibbsSampler(Dna, k, t, N)
def GibbsSampler(Dna, k, t, N):
    #randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
    Motifs = RandomMotifs(Dna, k, t)
    #BestMotifs ← Motifs
    BestMotifs = Motifs
    for i in range(1, N):
        #randomly generated integer between 1 and t
        r = random.randint(1, t)
        #Profile ← profile matrix formed from all strings in Motifs except for Motifi
        #Motifi ← Profile-randomly generated k-mer in the i-th string
        Motifi = Motifs[r-1]
        Motifs.remove(Motifs[r-1])
        Profile = ProfileWithPseudocounts(*Motifs)
        Motifs.insert(r-1, ProfileGeneratedString(Motifi, Profile, k))
        if Score(*Motifs) < Score(*BestMotifs):
                BestMotifs = Motifs
    return BestMotifs

