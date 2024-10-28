from unittest import TestCase
from alignment_algorithms import needleman_wunsch, smith_waterman
import numpy as np
char_to_index = {'A': 0, 'G': 1, 'C': 2, 'T': 3}

class TestNeedlemanWunsch(TestCase):

    def setUp(self):
        self.substitution_matrix = np.array([
            [1, -1, -1, -1],  # A
            [-1, 1, -1, -1],  # T
            [-1, -1, 1, -1],  # G
            [-1, -1, -1, 1]   # C
        ])
        self.char_to_index = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
        self.gap_penalty = -1

    def test_match(self):
        score_matrix = needleman_wunsch("A", "A", self.substitution_matrix, self.gap_penalty)
        self.assertEqual(score_matrix[1][1], 1)

    def test_mismatch(self):
        score_matrix = needleman_wunsch("A", "T", self.substitution_matrix, self.gap_penalty)
        self.assertEqual(score_matrix[1][1], -1)


    def test_alignment_with_gaps(self):
        score_matrix = needleman_wunsch("AGCT", "AGT", self.substitution_matrix, self.gap_penalty)
        self.assertEqual(score_matrix[4][3], 2)  # A, G match; C has a gap; T matches

    def test_longer_sequences(self):
        score_matrix = needleman_wunsch("AGCTG", "AGTG", self.substitution_matrix, self.gap_penalty)
        self.assertEqual(score_matrix[5][4], 3)  # Expected score based on the scoring matrix

    def test_longer_sequences_1(self):
        score_matrix = needleman_wunsch("AAAAGGGG", "AAAAGGGG", self.substitution_matrix, self.gap_penalty)
        self.assertEqual(score_matrix[8][8], 8)  # Expected score based on the scoring matrix

    def test_longer_sequences_2(self):
        score_matrix = needleman_wunsch("AAGGTTCC", "AAGGTCC", self.substitution_matrix, self.gap_penalty)
        self.assertEqual(score_matrix[8][7], 6)  # Expected score based on the scoring matrix


class TestSmithWaterman(TestCase):
    def setUp(self):
        self.substitution_matrix = np.array([
            [1, -1, -1, -1],
            [-1, 1, -1, -1],
            [-1, -1, 1, -1],
            [-1, -1, -1, 1]
        ])
        self.char_to_index = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
        self.gap_penalty = -1

    def test_match(self):
        score_matrix, max_positions, max_score = smith_waterman("A", "A", self.substitution_matrix, self.gap_penalty)
        self.assertEqual(max_score, 1)
        self.assertEqual(max_positions[0], (1, 1))

    def test_mismatch(self):
        score_matrix, max_positions, max_score = smith_waterman("A", "T", self.substitution_matrix, self.gap_penalty)
        self.assertEqual(max_score, 0)
        self.assertEqual(max_positions, [(1,1)])

    def test_local_alignment_with_gaps(self):
        score_matrix, max_positions, max_score = smith_waterman("AGTC", "AGT", self.substitution_matrix, self.gap_penalty)
        self.assertEqual(max_score, 3)
        self.assertIn((3, 3), max_positions)

    def test_empty_sequences(self):
        score_matrix, max_positions, max_score = smith_waterman("", "", self.substitution_matrix, self.gap_penalty)
        self.assertEqual(max_score, 0)
        self.assertEqual(max_positions, [])

    def test_one_empty_sequence(self):
        score_matrix, max_positions, max_score = smith_waterman("A", "", self.substitution_matrix, self.gap_penalty)
        self.assertEqual(max_score, 0)
        self.assertEqual(max_positions, [])

        score_matrix, max_positions, max_score = smith_waterman("", "A", self.substitution_matrix, self.gap_penalty)
        self.assertEqual(max_score, 0)
        self.assertEqual(max_positions, [])

    def test_identical_sequences(self):
        score_matrix, max_positions, max_score = smith_waterman("ATGC", "ATGC", self.substitution_matrix, self.gap_penalty)
        self.assertEqual(max_score, 4)
        self.assertIn((4, 4), max_positions)


