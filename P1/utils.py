from alignment_algorithms import needleman_wunsch, smith_waterman, backtrack


def validate_sequences(seq1: str, seq2: str, valid_chars: set = {'A', 'G', 'C', 'T'}) -> None:
    """
    Validate that the input sequences contain only allowed characters.

    Args:
        seq1: First DNA sequence (string).
        seq2: Second DNA sequence (string).
        valid_chars: A set of valid characters (e.g., {'A', 'G', 'C', 'T'}).

    Raises:
        ValueError: If any character in the sequences is not valid.
    """
    for seq in [seq1, seq2]:
        for char in seq:
            if char not in valid_chars:
                raise ValueError(f"Invalid character '{char}' found in sequence. Only 'A', 'G', 'C', 'T' are allowed.")


def find_global_alignments(seq1, seq2, substitution_matrix, gap_penalty, N):
    """
    Args:
        seq1: First DNA sequence (accepts only string made of 'A', 'T', 'G', 'C')
        seq2: Second DNA sequence (accepts only string made of 'A', 'T', 'G', 'C')
        N: Number of maximal optimal alignments
        substitution_matrix: A matrix that scores match/mismatch of two DNA sequences
        gap_penalty: A penalty of aligning '-' and non-empty element of the sequence

    Returns:
    alignments[:N]: best alignments
    global_score: corresponding score
    """
    validate_sequences(seq1, seq2)
    score_matrix = needleman_wunsch(seq1, seq2, substitution_matrix, gap_penalty)
    alignments = []
    backtrack(score_matrix, seq1, seq2, len(seq1), len(seq2), "", "", alignments, N, substitution_matrix, gap_penalty, True)
    global_score = score_matrix[len(seq1)][len(seq2)]
    return alignments[:N], global_score


def find_local_alignments(seq1, seq2, substitution_matrix, gap_penalty, N):
    """
    Args:
        seq1: First DNA sequence (accepts only string made of 'A', 'T', 'G', 'C')
        seq2: Second DNA sequence (accepts only string made of 'A', 'T', 'G', 'C')
        N: Number of maximal optimal alignments
        substitution_matrix: A matrix that scores match/mismatch of two DNA sequences
        gap_penalty: A penalty of aligning '-' and non-empty element of the sequence

    Returns:
    alignments[:N]: best alignments
    global_score: corresponding score
    """
    validate_sequences(seq1, seq2)
    score_matrix, max_positions, max_score = smith_waterman(seq1, seq2, substitution_matrix, gap_penalty)
    alignments = []
    for position in max_positions:
        backtrack(score_matrix, seq1, seq2, position[0], position[1], "", "", alignments, N, substitution_matrix, gap_penalty, False)
    return alignments[:N], max_score
