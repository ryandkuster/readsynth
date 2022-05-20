import os
import pandas as pd
import re
import sys
import unittest

sys.path.append("../")
import readsynth as rs
import scripts.digest_genomes as dg
import scripts.digest_genomes_iso as dgi
import scripts.prob_n_copies as pnc
import scripts.prob_n_copies_iso as pnci

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

    def test_open_fastq_1(self):
        pass

    def test_check_genomes_1(self):
        '''
        given an input df of genome file locations and abundances for
        each, return df of file locations with abundances identical to
        the expected file saved as example_metagenome_df.csv
        '''
        args = Variables()
        args.genome = 'test_data/example_metagenome.csv'
        expected_df = pd.read_csv('test_data/example_metagenome_df.csv')
        df = rs.check_genomes(args.genome)
        self.assertTrue(expected_df.equals(df))

    def test_check_genomes_2(self):
        '''
        given an input df of genome file locations and abundances for
        each, return df of file locations with abundances identical to
        the expected file saved as example_metagenome_df.csv
        '''
        args = Variables()
        args.genome = 'test_project3/metagenome_file.csv'
        expected_df = pd.read_csv('test_project3/metagenome_file_df.csv')
        df = rs.check_genomes(args.genome)
        self.assertTrue(expected_df.equals(df))


class TestIsoFunctions(unittest.TestCase):

    def test_iupac_motifs_iso_1(self):
        '''
        given 'NN/NNNNNNNNNNCGANNNNNNTGCNNNNNNNNNNNN/'
        expected:
            dictionary with sequence as key
            index of '/' is the value
        '''
        args = Variables()
        args.iso = ['NN/NNNNNNNNNNCGANNNNNNTGCNNNNNNNNNNNN/']
        expected = {'[ACGT][ACGT][ACGT][ACGT][ACGT][ACGT][ACGT][ACGT][ACGT]' +
                    '[ACGT][ACGT][ACGT]CGA[ACGT][ACGT][ACGT][ACGT][ACGT]' +
                    '[ACGT]TGC[ACGT][ACGT][ACGT][ACGT][ACGT][ACGT][ACGT]' +
                    '[ACGT][ACGT][ACGT][ACGT][ACGT]': [2, 36],
                    '[ACGT][ACGT][ACGT][ACGT][ACGT][ACGT][ACGT][ACGT][ACGT]' +
                    '[ACGT][ACGT][ACGT]GCA[ACGT][ACGT][ACGT][ACGT][ACGT]' +
                    '[ACGT]TCG[ACGT][ACGT][ACGT][ACGT][ACGT][ACGT][ACGT]' +
                    '[ACGT][ACGT][ACGT][ACGT][ACGT]': [0, 34]}
        motif_dt = rs.iupac_motifs_iso(args.iso)
        self.assertEqual(motif_dt, expected)

    def test_dgi_digest_seq_1(self):
        args = Variables()
        seq = 'A'*22 + 'CGAAAAAAATGC' + 'A'*66
        motif_dt = {'[ACGT]'*12 + 'CGA' + '[ACGT]'*6 + 'TGC' + '[ACGT]'*12:\
                    [2, 36],
                    '[ACGT]'*12 + 'CGA' + '[ACGT]'*6 + 'TGC' + '[ACGT]'*12:\
                    [0, 34]}
        expected = [['AAAAAAAAAAAACGAAAAAAATGCAAAAAAAAAAAA', 10, 46,
                     '[ACGT]'*12 + 'CGA' + '[ACGT]'*6 + 'TGC' + '[ACGT]'*12,
                     '[ACGT]'*12 + 'CGA' + '[ACGT]'*6 + 'TGC' + '[ACGT]'*12,
                     0]]
        seq_ls = dgi.digest_seq(0, seq, motif_dt, 200)
        self.assertEqual(seq_ls, expected)

    def test_dgi_get_query_len_1(self):
        args = Variables()
        query = '[ACGT]'*12 + 'CGA' + '[ACGT]'*6 + 'TGC' + '[ACGT]'*12
        mot_len = dgi.get_query_len(query)
        self.assertEqual(mot_len, 36)

    def test_pnci_find_overlaps_1(self):
        args = Variables()
        df = pd.read_csv('test_data/iso_genomes/raw_digest_test1.fasta.csv')
        df = df.rename(columns={'internal': 'overlaps'})
        expected = [None,
                    [2], [1],
                    [4, 5], [3, 5], [3, 4],
                    None,
                    [8], [7, 9], [8],
                    [11, 12], [10, 12], [10, 11, 13], [12],
                    [15], [14, 16], [15, 17], [16, 18], [17]]
        df = pnci.find_overlaps(df)
        self.assertEqual(df['overlaps'].tolist(), expected)

    def test_pnci_count_events_1(self):
        linked = [[2], [1]]
        expected = [(1,), (2,)]
        first, last, events = pnci.count_events(linked)
        self.assertEqual((first, last), (1, 2))
        self.assertCountEqual(events, expected)

    def test_pnci_count_events_2(self):
        linked = [[4, 5], [3, 5], [3, 4]]
        expected = [(3,), (4,), (5,)]
        first, last, events = pnci.count_events(linked)
        self.assertEqual((first, last), (3, 5))
        self.assertCountEqual(events, expected)

    def test_pnci_count_events_3(self):
        linked = [[8],[7,9],[8]]
        expected = [(7, 9), (8,)]
        first, last, events = pnci.count_events(linked)
        self.assertEqual((first, last), (7, 9))
        self.assertCountEqual(events, expected)

    def test_pnci_count_events_4(self):
        linked = [[11, 12], [10, 12], [10, 11, 13], [12]]
        expected = [(10, 13), (11, 13), (12,)]
        first, last, events = pnci.count_events(linked)
        self.assertEqual((first, last), (10, 13))
        self.assertCountEqual(events, expected)

    def test_pnci_count_events_5(self):
        linked = [[15], [14, 16], [15, 17], [16, 18], [17]]
        expected = [(14, 16, 18), (14, 17), (15,17), (15, 18)]
        first, last, events = pnci.count_events(linked)
        self.assertEqual((first, last), (14, 18))
        self.assertCountEqual(events, expected)

    def test_pnci_apply_probs_1(self):
        events = [(3,), (4,), (5,)]
        expected = [1/3, 1/3, 1/3]
        cluster_probs = pnci.apply_probs(3, 5, events)
        self.assertEqual(cluster_probs, expected)

    def test_pnci_apply_probs_2(self):
        '''
        4 possible events: 14 in 2, 15 in 2, 16 in 1, 17 in 2, 18 in 2
        '''
        events = [(14, 16, 18), (14, 17), (15,17), (15, 18)]
        expected = [2/4, 2/4, 1/4, 2/4, 2/4]
        cluster_probs = pnci.apply_probs(14, 18, events)
        self.assertEqual(cluster_probs, expected)

    def test_pnci_calculate_prob(self):
        df = pd.read_csv('test_data/iso_genomes/raw_digest_test1.fasta.csv')
        df = df.rename(columns={'internal': 'overlaps'})
        df = pnci.find_overlaps(df)
        df = pnci.calculate_prob(df)
        expected = [1,
                    1/2, 1/2,
                    1/3, 1/3, 1/3,
                    1,
                    1/2, 1/2, 1/2,
                    1/3, 1/3, 1/3, 2/3,
                    1/2, 1/2, 1/4, 1/2, 1/2]
        self.assertEqual(df['probability'].tolist(), expected)

    def test_pnci_get_len_freqs(self):
        df = pd.read_csv('test_data/iso_genomes/counts_test1.fasta.csv')
        expected = {i:0 for i in range(0, 34)}
        expected[34] = sum(df['adj_prob'].tolist())
        df = pnci.get_len_freqs(df, 34)
        self.assertEqual(df, expected)


class TestProcessing(unittest.TestCase):
    '''
    smoke tests for simple processing pipeline examples

    tests to perform:
        identical RE for m1/m2
        different RE for m1/m2
        nested RE for m1/m2 (identical)
        nested RE for m1/m2 (different)
        iso RE cutter (standard cut regions)
        iso RE cutter (overlapping cut regions)
    '''

#    def test_process_genomes_1(self):
#        '''
#        '''
#        args = Variables()
#        args.genome = 'test_project3/metagenome_file.csv'
#        args.o = 'test_project3'
#        args.motif_dt = {'TTAA': 1, 'GCGC': 3}
#        args.motif_dt1 = {'TTAA': 1}
#        args.motif_dt2 = {'GCGC': 3}
#        args.cut_prob = 1
#        genomes_df = pd.read_csv('test_project3/metagenome_file_df.csv')
#        genomes_df, total_freqs = rs.process_genomes(args, genomes_df)
##        self.assertEqual(args.m1[0], expected)

#    def test_process_genomes_2(self):
#        '''
#        '''
#        args = Variables()
#        args.genome = 'test_project4/metagenome_file.csv'
#        args.o = 'test_project4'
#        args.motif_dt = {'GCGC': 3, 'TTAA': 1}
#        args.motif_dt1 = {'GCGC': 3}
#        args.motif_dt2 = {'TTAA': 1}
#        args.cut_prob = 0.5
#        genomes_df = pd.read_csv('test_project4/metagenome_file_df.csv')
#        genomes_df, total_freqs = rs.process_genomes(args, genomes_df)
#        df = pd.read_csv(
#            'test_project4/raw_digest_Escherichia_coli_plasmid.fasta.csv')
##        self.assertEqual(args.m1[0], expected)


    def test_reverse_comp_1(self):
        '''
        test basic reverse complementarity function
        '''
        seq = 'ACGTN'
        expected = 'NACGT'
        new = rs.reverse_comp(seq)
        self.assertEqual(new, expected)

    def test_reverse_comp_2(self):
        '''
        test redundant IUPAC reverse complementarity function
        '''
        seq = 'RVDAB'
        expected = 'VTHBY'
        new = rs.reverse_comp(seq)
        self.assertEqual(new, expected)

    def test_reverse_comp_3(self):
        '''
        test highly redundant IUPAC reverse complementarity function
        '''
        seq = 'ACYGRYSCCTATTACNNCTNVBBVAYSR'
        expected = 'YSRTBVVBNAGNNGTAATAGGSRYCRGT'
        new = rs.reverse_comp(seq)
        self.assertEqual(new, expected)


class TestDigestGenomes(unittest.TestCase):

    def test_dg_main_1(self):
        '''
        digest test_project1/fake_genome_1.fasta and confirm number of
        internal cut sites behaves as expected
        '''
        args = Variables()
        args.genome = 'test_project1/fake_genome_1.fasta'
        args.motif_dt = {'TTAA': 1}
        df = dg.main(args)
        self.assertEqual(df.iloc[df[(df['start'] == 25) & (df['end'] == 125)]\
            .index]['internal'].to_list()[0], 3)

    def test_dg_main_2(self):
        pass
        args = Variables()
        args.genome = 'test_project2/fake_genome_2.fasta'
        args.motif_dt = {'TTAA': 1, 'GCGC': 3}
        df = dg.main(args)
        self.assertEqual(df.iloc[df[(df['start'] == 25) & (df['end'] == 75)]\
            .index]['internal'].to_list()[0], 1)


if __name__ == "__main__":
    unittest.main()
