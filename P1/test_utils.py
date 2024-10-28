import unittest
from unittest import TestCase

from utils import validate_sequences

class TestValidateSequences(unittest.TestCase):
    def test_valid_sequences(self):
        try:
            validate_sequences("AGCT", "CGTA")
        except ValueError:
            self.fail("validate_sequences raised ValueError unexpectedly!")

    def test_invalid_character_in_first_sequence(self):
        with self.assertRaises(ValueError) as context:
            validate_sequences("AGXT", "CGTA")
        self.assertEqual(str(context.exception), "Invalid character 'X' found in sequence. Only 'A', 'G', 'C', 'T' are allowed.")

    def test_invalid_character_in_second_sequence(self):
        with self.assertRaises(ValueError) as context:
            validate_sequences("AGCT", "C1GT")
        self.assertEqual(str(context.exception), "Invalid character '1' found in sequence. Only 'A', 'G', 'C', 'T' are allowed.")

    def test_both_sequences_invalid(self):
        with self.assertRaises(ValueError) as context:
            validate_sequences("AGXT", "C1GT")
        self.assertEqual(str(context.exception), "Invalid character 'X' found in sequence. Only 'A', 'G', 'C', 'T' are allowed.")

    def test_empty_sequences(self):
        try:
            validate_sequences("", "")
        except ValueError:
            self.fail("validate_sequences raised ValueError unexpectedly for empty sequences!")

    def test_all_valid_chars(self):
        try:
            validate_sequences("AGCTAGC", "TCGATAG")
        except ValueError:
            self.fail("validate_sequences raised ValueError unexpectedly for valid sequences with multiple characters!")

    def test_invalid_sequence_with_special_chars(self):
        with self.assertRaises(ValueError) as context:
            validate_sequences("AGC$", "TGCA")
        self.assertEqual(str(context.exception), "Invalid character '$' found in sequence. Only 'A', 'G', 'C', 'T' are allowed.")
