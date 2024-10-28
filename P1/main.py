from utils import find_local_alignments, find_global_alignments
import pandas as pd
import sys

if __name__ == "__main__":
    # parameters
    FILEPATH = "data/substitution_matrix.csv"
    OUTPUT_FILEPATH = "data/output.log"
    substitution_matrix = pd.read_csv(FILEPATH).to_numpy()
    gap_penalty = -2
    n = 3

    # Example usage
    seq1 = "TATA"
    seq2 = "ATAT"

    original_stdout = sys.stdout
    with open(OUTPUT_FILEPATH, 'w') as f:
        sys.stdout = f

        # Find n maximal global alignments
        global_alignments, global_score = find_global_alignments(seq1, seq2, substitution_matrix, gap_penalty, n)
        for i, alignment in enumerate(global_alignments):
            print(f"Global alignment no {i+1}:\n{alignment[0]}\n{alignment[1]}\nScore: {global_score}\n")

        # Find n maximal local alignments
        local_alignments, local_score = find_local_alignments(seq1, seq2, substitution_matrix, gap_penalty, n)
        for i, alignment in enumerate(local_alignments):
            print(f"Local alignment no {i+1}:\n{alignment[0]}\n{alignment[1]}\nScore: {local_score}\n")



