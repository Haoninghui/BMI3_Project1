import math


k = 0.1
lam = 0.9  # range 0.8-0.9


def SW_score_filter(score, threshold):
    return score <= threshold


def Escore_filter(SW_score, m, n, Escore_threshold):
    constant = k * m * n
    E = constant * math.exp(-lam * SW_score)
    return E >= Escore_threshold
