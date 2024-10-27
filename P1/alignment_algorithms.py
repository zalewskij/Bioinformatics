import numpy as np

char_to_index = {'A': 0, 'G': 1, 'C': 2, 'T': 3}


def needleman_wunsch(seq1: str, seq2: str, substitution_matrix: np.array, gap_penalty: int):
    """
    Args:
        seq1: First DNA sequence (accepts only string made of 'A', 'T', 'G', 'C')
        seq2: Second DNA sequence (accepts only string made of 'A', 'T', 'G', 'C')
        substitution_matrix: A matrix that scores match/mismatch of two DNA sequences
        gap_penalty: A penalty of aligning '-' and non-empty element of the sequence

    Returns:
        score_matrix: Filled matrix of scores (to get th alignment need to backtrack)
    """

    len_seq1 = len(seq1)
    len_seq2 = len(seq2)

    # Initialize score matrix (aligning characters with gaps)
    score_matrix = np.zeros((len_seq1 + 1, len_seq2 + 1), dtype=int)
    for i in range(len_seq1 + 1):
        score_matrix[i][0] = gap_penalty * i
    for j in range(len_seq2 + 1):
        score_matrix[0][j] = gap_penalty * j

    # Compute score for subsequent types of modifications by filling the scoring matrix
    for i in range(1, len_seq1 + 1):
        for j in range(1, len_seq2 + 1):
            char1_index = char_to_index[seq1[i - 1]]
            char2_index = char_to_index[seq2[j - 1]]
            match = score_matrix[i - 1][j - 1] + substitution_matrix[char1_index][char2_index]
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty
            score_matrix[i][j] = max(match, delete, insert)

    return score_matrix


def smith_waterman(seq1, seq2, substitution_matrix, gap_penalty):
    """
        Args:
            seq1: First DNA sequence (accepts only string made of 'A', 'T', 'G', 'C')
            seq2: Second DNA sequence (accepts only string made of 'A', 'T', 'G', 'C')
            substitution_matrix: A matrix that scores match/mismatch of two DNA sequences
            gap_penalty: A penalty of aligning '-' and non-empty element of the sequence

        Returns:
            score_matrix: Filled matrix of scores (to get th alignment need to backtrack)
            max_positions: Positions on scoring matrix with maximal score
            max_score: Maximal score
        """
    len_seq1 = len(seq1)
    len_seq2 = len(seq2)

    # Initialize score matrix (no aligning characters with gaps)
    score_matrix = np.zeros((len_seq1 + 1, len_seq2 + 1), dtype=int)
    max_score = 0
    max_positions = []

    # Compute score for subsequent types of modifications by filling the scoring matrix
    for i in range(1, len_seq1 + 1):
        for j in range(1, len_seq2 + 1):
            char1_index = char_to_index[seq1[i - 1]]
            char2_index = char_to_index[seq2[j - 1]]

            match = score_matrix[i - 1][j - 1] + substitution_matrix[char1_index][char2_index]
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty

            # We can reset to 0 for local alignment (sequence can start at any point)
            score_matrix[i][j] = max(0, match, delete, insert)

            # Track the maximum score positions
            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_positions = [(i, j)]
            elif score_matrix[i][j] == max_score:
                max_positions.append((i, j))

    return score_matrix, max_positions, max_score


def backtrack(score_matrix: np.array, seq1: str, seq2: str, i: int, j: int, alignment1: str, alignment2: str,
              alignments: list, N: int, substitution_matrix: np.array, gap_penalty: int, global_alignment: bool = True):
    """
    Args:
        score_matrix: A score matrix from needleman_wunsh or smith_waterman algorithm
        seq1: First DNA sequence (accepts only string made of 'A', 'T', 'G', 'C')
        seq2: Second DNA sequence (accepts only string made of 'A', 'T', 'G', 'C')
        i: Current index in the scoring matrix (horizontal)
        j: Current index in the scoring matrix (vertical)
        alignment1: modified seq1
        alignment2: modified seq2
        alignments: A list of up to N collected alignments
        N: Number of maximal optimal alignments
        substitution_matrix: A matrix that scores match/mismatch of two DNA sequences
        gap_penalty: A penalty of aligning '-' and non-empty element of the sequence
        global_alignment: A boolean variable that indicates if we look for global or local alignment

    Returns:
    """

    # base: For global and local alignments - stop if already found N maximal alignments
    if len(alignments) >= N:
        return

    # base: For local alignment - stop if the score is zero
    if not global_alignment and score_matrix[i][j] == 0:
        alignments.append((alignment1[::-1], alignment2[::-1]))
        return

    # base: for global alignment - stop if we've reached the start (0, 0)
    if global_alignment and i == 0 and j == 0:
        alignments.append((alignment1[::-1], alignment2[::-1]))
        return

    # Character indices for current position
    char1_index = char_to_index[seq1[i - 1]] if i > 0 else None
    char2_index = char_to_index[seq2[j - 1]] if j > 0 else None

    # Matching
    if i > 0 and j > 0:
        match_score = score_matrix[i - 1][j - 1] + substitution_matrix[char1_index][char2_index]
        if score_matrix[i][j] == match_score:
            backtrack(score_matrix, seq1, seq2, i - 1, j - 1,
                      alignment1 + seq1[i - 1], alignment2 + seq2[j - 1],
                      alignments, N, substitution_matrix, gap_penalty, global_alignment)

    # Delete
    if i > 0:
        delete_score = score_matrix[i - 1][j] + gap_penalty
        if score_matrix[i][j] == delete_score:
            backtrack(score_matrix, seq1, seq2, i - 1, j,
                      alignment1 + seq1[i - 1], alignment2 + "-",
                      alignments, N, substitution_matrix, gap_penalty, global_alignment)

    # Insert
    if j > 0:
        insert_score = score_matrix[i][j - 1] + gap_penalty
        if score_matrix[i][j] == insert_score:
            backtrack(score_matrix, seq1, seq2, i, j - 1,
                      alignment1 + "-", alignment2 + seq2[j - 1],
                      alignments, N, substitution_matrix, gap_penalty, global_alignment)
