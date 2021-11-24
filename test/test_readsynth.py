import os
import pandas as pd
import re
import sys
import unittest

sys.path.append("../")
import readsynth as rs
import digest_genomes as dg

'''
usage:
python3 -m unittest test_readsynth.py
'''

class Variables():

    def __init__(self):
        self.genome = 'filler'
        self.m1 = ['T/TAA']
        self.m2 = ['T/TAA']
        self.l = 150
        #self.test
        #self.t
        #self.n
        self.mean = 400
        self.up_bound = 500
        self.min = 6
        self.max = 700
        #self.cut_prob
        #self.a1
        #self.a2
        #self.a1s
        #self.a2s
        #self.q1
        #self.q2
        #self.r1
        #self.r2
        #self.p


class TestInitFunctions(unittest.TestCase):

    def test_check_for_enzymes_1(self):
        '''
        open pickled dictionary of enzymes
        return 'A/CCGGT' for 'AGEI'
        '''
        args = Variables()
        args.m1 = ['AGEI']
        expected = 'A/CCGGT'
        args.m1, args.m2 = rs.check_for_enzymes(args)
        self.assertEqual(args.m1[0], expected)

    def test_check_for_enzymes_2(self):
        '''
        open pickled dictionary of enzymes
        return 'A/CCGGT' for 'agei'
        '''
        args = Variables()
        args.m1 = ['agei']
        expected = 'A/CCGGT'
        args.m1, args.m2 = rs.check_for_enzymes(args)
        self.assertEqual(args.m1[0], expected)

    def test_check_for_enzymes_3(self):
        '''
        open pickled dictionary of enzymes
        return 'A/CCGGT' for 'A/CCGGT'
        '''
        args = Variables()
        args.m1 = ['A/CCGGT']
        expected = 'A/CCGGT'
        args.m1, args.m2 = rs.check_for_enzymes(args)
        self.assertEqual(args.m1[0], expected)

    def test_iupac_motifs_1(self):
        '''
        given 'A/CCGGT'
        expected:
            dictionary with sequence as key
            index of '/' is the value
        '''
        args = Variables()
        args.m1 = ['A/CCGGT']
        expected = {'ACCGGT': 1}
        motif_dt = rs.iupac_motifs(args.m1)
        self.assertEqual(motif_dt, expected)

    def test_iupac_motifs_2(self):
        '''
        given 'C/YCGRG'
        expected:
            dictionary with sequence as key
            index of '/' is the value
            Y and R redundancy captured in regex form
        '''
        args = Variables()
        args.m1 = ['C/YCGRG']
        expected = {'C[CT]CG[AG]G': 1}
        motif_dt = rs.iupac_motifs(args.m1)
        self.assertEqual(motif_dt, expected)

    def test_get_adapters_1(self):
        pass

    def test_get_sbs_start_1(self):
        pass

    def test_get_qscores_1(self):
        pass

    def test_open_fastq_1(self):
        pass

    def test_check_genomes_1(self):
        '''
        given an input df of genome file locations and abundances for
        each, return df of file locations with abundances identical to
        the expected file saved as example_metagenome_df.csv
        '''
        args = Variables()
        args.genome = 'test_data/example_metagenome_df.csv'
        expected_df = pd.read_csv(args.genome)
        df = rs.check_genomes('test_data/example_metagenome.csv')
        self.assertTrue(expected_df.equals(df))


class TestProcessing(unittest.TestCase):
#
#    def test_process_genomes_1(self):
#        '''
#        '''
#        args = Variables()
#        args.genome = 
#        args.proj = 
#        genomes_df = pd.DataFrame(np.array([['a',1]]))
#        genomes_df, total_freqs = rs.process_genomes(args, genomes_df)
#        self.assertEqual(args.m1[0], expected)

    def test_reverse_comp_1(self):
        '''
        '''
        seq = 'ACGTN'
        expected = 'NACGT'
        new = rs.reverse_comp(seq)
        self.assertEqual(new, expected)

    def test_reverse_comp_2(self):
        '''
        '''
        seq = 'RVDAB'
        expected = 'VTHBY'
        new = rs.reverse_comp(seq)
        self.assertEqual(new, expected)

    def test_reverse_comp_3(self):
        '''
        '''
        seq = 'ACYGRYSCCTATTACNNCTNVBBVAYSR'
        expected = 'YSRTBVVBNAGNNGTAATAGGSRYCRGT'
        new = rs.reverse_comp(seq)
        self.assertEqual(new, expected)


class TestDigestGenomes(unittest.TestCase):

    def test_dg_main_1(self):
        pass
        args = Variables()
        args.genome = 'test_project1/fake_genome_1.fasta'
        motif_dt = {'TTAA': 1}
        df = dg.main(motif_dt, args)
        self.assertEqual(df.iloc[df[(df['start'] == 25) & (df['end'] == 125)]\
                         .index]['internal'].to_list()[0], 3)

    def test_dg_main_2(self):
        pass
        args = Variables()
        args.genome = 'test_project2/fake_genome_2.fasta'
        motif_dt = {'TTAA': 1, 'GCGC': 3}
        df = dg.main(motif_dt, args)









if __name__ == "__main__":
    unittest.main()
