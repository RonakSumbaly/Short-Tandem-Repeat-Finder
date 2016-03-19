import re
import improved_helpers
import logging as logger

logger.basicConfig(level=logger.INFO, format='> %(message)s')

match_award = 10
mismatch_penalty = -5
gap_penalty = -5


def make_matrix(shape):
    return [[0] * shape[1] for i in xrange(shape[0])]


def match_score(alpha, beta):
    if alpha == beta:
        return match_award
    elif alpha == '-' or beta == '-':
        return gap_penalty
    else:
        return mismatch_penalty


def add_symbols(align_seq_1, align_seq_2):
    align_seq_1 = align_seq_1[::-1]  # reverse sequence 1
    align_seq_2 = align_seq_2[::-1]  # reverse sequence 2

    # calculate identity, score and aligned sequeces
    symbol = ''
    score = 0
    identity = 0
    for i in range(0, len(align_seq_1)):
        # if two AAs are the same, then output the letter
        if align_seq_1[i] == align_seq_2[i]:
            symbol = symbol + align_seq_1[i]
            identity += 1
            score += match_score(align_seq_1[i], align_seq_2[i])

        # if they are not identical and none of them is gap
        elif align_seq_1[i] != align_seq_2[i] and align_seq_1[i] != '-' and align_seq_2[i] != '-':
            score += match_score(align_seq_1[i], align_seq_2[i])
            symbol += ' '

        # if one of them is a gap, output a space
        elif align_seq_1[i] == '-' or align_seq_2[i] == '-':
            symbol += ' '
            score += gap_penalty

    identity = float(identity) / len(align_seq_1) * 100

    return align_seq_1, align_seq_2


def waterman_algorithm(seq_1, seq_2):
    global max_i, max_j

    m, n = len(seq_1), len(seq_2)  # length of two sequences

    # Generate DP table and traceback path pointer matrix
    score = make_matrix((m + 1, n + 1))  # the DP table
    pointer = make_matrix((m + 1, n + 1))  # to store the traceback path

    max_score = 0  # initial maximum score in DP table

    # Calculate DP table and mark pointers
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            score_diagonal = score[i - 1][j - 1] + match_score(seq_1[i - 1], seq_2[j - 1])
            score_up = score[i][j - 1] + gap_penalty
            score_left = score[i - 1][j] + gap_penalty
            score[i][j] = max(0, score_left, score_up, score_diagonal)
            if score[i][j] == 0:
                pointer[i][j] = 0  # 0 means end of the path
            if score[i][j] == score_left:
                pointer[i][j] = 1  # 1 means trace up
            if score[i][j] == score_up:
                pointer[i][j] = 2  # 2 means trace left
            if score[i][j] == score_diagonal:
                pointer[i][j] = 3  # 3 means trace diagonal
            if score[i][j] >= max_score:
                max_i = i
                max_j = j
                max_score = score[i][j]

    align_1, align_2 = '', ''  # initial sequences

    i, j = max_i, max_j  # indices of path starting point

    # traceback, follow pointers
    try:
        while pointer[i][j] != 0:
            if pointer[i][j] == 3:
                align_1 += seq_1[i - 1]
                align_2 += seq_2[j - 1]
                i -= 1
                j -= 1
            elif pointer[i][j] == 2:
                align_1 += '-'
                align_2 += seq_2[j - 1]
                j -= 1
            elif pointer[i][j] == 1:
                align_1 += seq_1[i - 1]
                align_2 += '-'
                i -= 1
    except:
        return False, False

    return add_symbols(align_1, align_2)


def check_for_insertions(align_seq_1, align_seq_2, start_ref, check_ref):
    """
    :return: deletions (function_name incorrect - cannot change as will have to change in lots of places :P)
    """
    indices = [(i.start(), i.end()) for i in re.finditer('-+', align_seq_1)]

    if len(indices) == 1:
        check = align_seq_2[0:indices[0][0]]
        i = check_ref.find(check)
        start_pos = (i + len(check) + start_ref)
        DEL = align_seq_2[indices[0][0]:indices[0][1]]
        improved_helpers.add_deletions(start_pos, DEL)
    elif len(indices) == 2 and abs(indices[0][0] - indices[1][0]) < 6:
        if abs(indices[0][0] - indices[0][1]) == 1 and abs(indices[1][0] - indices[1][1]) == 1:
            check = align_seq_2[0:indices[0][0]]
            i = check_ref.find(check)
            DEL = align_seq_2[indices[0][0] + 1:indices[1][0] + 1]
            start_pos = (i + len(check) + start_ref + len(DEL) - 1)
            improved_helpers.add_deletions(start_pos, DEL)  # made a mistake but


def check_for_deletions(align_seq_1, align_seq_2, start_ref, check_ref):
    """
    :return: insertions (function_name incorrect - cannot change as will have to change in lots of places :P)
    """
    indices = [(i.start(), i.end()) for i in re.finditer('-+', align_seq_2)]

    if len(indices) == 1:
        check = align_seq_1[0:indices[0][0]]
        i = check_ref.find(check)
        start_pos = (i + len(check) + start_ref)
        INS = align_seq_1[indices[0][0]:indices[0][1]]
        improved_helpers.add_insertions(start_pos, INS)

    elif len(indices) == 2 and abs(indices[0][0] - indices[1][0]) < 6:
        if abs(indices[0][0] - indices[0][1]) == 1 and abs(indices[1][0] - indices[1][1]) == 1:
            check = align_seq_1[0:indices[0][0]]
            i = check_ref.find(check)
            INS = align_seq_1[indices[0][0] + 1:indices[1][0] + 1]
            start_pos = (i + len(check) + start_ref + len(INS) - 1)
            improved_helpers.add_insertions(start_pos, INS)