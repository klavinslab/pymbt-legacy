from collections import Counter


def repeat_check(seq, n):
    n_mers = [seq[i:i + n] for i in range(len(seq) - n)]
    counted = Counter(n_mers)
    return counted.most_common()