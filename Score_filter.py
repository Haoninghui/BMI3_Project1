import math


k = 0.1
lam = 0.9  # range 0.8-0.9


def SW_score_filter(score: int, threshold: int) -> bool:
    """
    Judge if SW score is below threshold.
    :param score: SW score calculated.
    :param threshold: Set threshold.
    :return: Below threshold or above threshold.
    """
    return score <= threshold


def Escore_filter(SW_score: int, m: int, n: int, Escore_threshold: float) -> bool:
    """
    Judge if E-score is above threshold
    :param SW_score: Smith-Waterman score.
    :param m: Query length.
    :param n: Reference length.
    :param Escore_threshold: Set threshold for E-score.
    :return: Above threshold or below threshold.
    """
    constant = k * m * n
    E = constant * math.exp(-lam * SW_score)
    return E >= Escore_threshold
