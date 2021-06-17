import unittest
import re
import sys
import os
sys.path.append("../")
#from readsynth import orientation_test


'''
usage:
python3 -m unittest test_readsynth.py
'''

class TestReadsynth(unittest.TestCase):

    def test_orientation_test_1(self):
        '''
        test regex function
        '''
        seq = 'CATGAAAAAAAAAACATG'
        motif_dt1 = {'CA[AT]G': 4}
        motif_dt2 = {'CA[AT]G': 4}
        expected = 'AAAAAAAAAA'
        result = orientation_test(seq, motif_dt1, motif_dt2)
        self.assertEqual(result, expected)

    def test_orientation_test_2(self):
        '''
        test regex function
        '''
        seq = 'CAAGAAAAAAAAAACATG'
        motif_dt1 = {'CA[AT]G': 4}
        motif_dt2 = {'CA[AT]G': 4}
        expected = 'AAAAAAAAAA'
        result = orientation_test(seq, motif_dt1, motif_dt2)
        self.assertEqual(result, expected)

    def test_orientation_test_3(self):
        '''
        test regex function
        '''
        seq = 'CACGAAAAAAAAAACATG'
        motif_dt1 = {'CA[AT]G': 4}
        motif_dt2 = {'CA[AT]G': 4}
        expected = 'AAAAAAAAAA'
        result = orientation_test(seq, motif_dt1, motif_dt2)
        self.assertEqual(result, None)


def orientation_test(seq, motif_dt1, motif_dt2):
    '''
    this is a hacky way to test this function
    but the original code uses several global variables
    that are not easily mocked in unittest
    '''
    for motif1, offset1 in motif_dt1.items():
        if re.search('^'+motif1, seq):
            for motif2, offset2 in motif_dt2.items():
                if re.search(motif2+'$', seq):
                    return seq[offset1:-offset2]

    return None


if __name__ == "__main__":
    unittest.main()
