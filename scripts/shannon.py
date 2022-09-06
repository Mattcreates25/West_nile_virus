import collections

from scipy.stats import entropy

def estimate_shannon_entropy(dna_sequence):
    bases = collections.Counter([tmp_base for tmp_base in dna_sequence])
    # define distribution
    dist = [x/sum(bases.values()) for x in bases.values()]

    # use scipy to calculate entropy
    entropy_value = entropy(dist, base=2)

    return entropy_value


print(estimate_shannon_entropy('results/all_alignment.fasta'))

