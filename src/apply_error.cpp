/*
 * File: apply_error.cpp
 * Author: Ryan Kuster
 * Date: March 4, 2024
 * Description:
 *   This script applies SNP-only, Phred-like error to input fastq-formatted
 *   sequence reads.
 * License: Apache-2.0 license
 */

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <vector>


char mutate_base (char base) {
    /*
     * create "bases" string with anticipated nucleotides
     * remove the current base to be mutated from "bases"
     * create a string array with the 3 non-input bases
     * randomly choose a non-input base
     * return the mutated base (mut_base)
     */
    std::string bases = "ACGT";
    bases.erase(std::remove(bases.begin(), bases.end(), base), bases.end());
    std::string choices[3];

    for (int i = 0; i < 3; i++) {
        choices[i] = bases[i];
    }

    std::string tmp_mut = choices[rand() % 3];
    assert(tmp_mut.size() == 1);
    char mut_base = tmp_mut[0];
    return mut_base;
}


bool miscall_test(double prob){
    /*
     * return true if rand() less than the current q score probability
     * multiplied by the maximum random constant (+ 1 for inclusivity)
     *
     * example:
     * if prob = 1, RAND_MAX + 1 ensures rand() cannot be greater so
     * false will always be returned (as expected)
     */
    return rand() < prob * (RAND_MAX+1.0);
}


std::string mutate_seq (std::string seq, std::string score, double scores_arr[], std::string phred_scores) {
    /*
     * initialize empty output string for mutated sequence (mut_seq)
     * iterate over each base in input (seq)
     */
    std::string mut_seq = "";
    for (int i = 0; i < seq.length(); i++) {
        char base = seq[i];
        
        /*
         * if base is N, append to mut_seq and continue
         */
        if (base == 'N') {
            mut_seq = mut_seq + base;
            continue;
        }
        
        /*
         * index the phred_scores array for the current q score position
         * look up the probability of the q score in scores_arr
         * miscall returns bool, determining if current base (seq[i])
         * will be mutated using mutate_base function
         */
        std::size_t idx = phred_scores.find(score[i]);
        double prob = scores_arr[idx];
        bool miscall = miscall_test (prob);

        if (miscall == true) {
            char mut_base = mutate_base (base);
            mut_seq = mut_seq + mut_base;
        }
        else {
            mut_seq = mut_seq + base;
        }
    }

    return mut_seq;
}


int main (int argc, char **argv)
{
    /*
     * apply_error takes an existing fastq-formatted file and performs
     * per-base (SNP only) error simulation applying the q score
     * probability at face-value
     *
     * input fastq file (file_name):
     *     contain reference genome sequences with no polymorphisms
     * output fastq file name (out_name)
     *     output fastq file to write
     */

    std::string file_name = argv[1];
    std::string out_name = argv[2];

    std::ifstream file(file_name);
    std::ofstream out(out_name);

    /*
     * create a string of possible phred score ascii characters
     * create an empty array of length N (score_arr)
     * populate scores_arr with probabilities ( 10^(-q/10) )
     */

    std::string phred_scores = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ";
    int N = phred_scores.length();
    double scores_arr[N];

    for (int i = 0; i < N; i++) {
        double y = -(double)i/10;
        scores_arr[i] = pow(10, y);
    }

    /*
     * initialize strings for each fastq read element
     */

    int line_no = 0;
    std::string line;
    std::string header = "";
    std::string seq = "";
    std::string plus = "+";
    std::string score = "";

    srand(time(0));

    /*
     * iterate through input fastq file
     * store each element (header, seq, score)
     * at each fourth line, mutate the sequence using the corresponding
     * q scores from each line (mutate_seq function)
     * write the four fastq elements to the "out" file
     */
    while (std::getline(file, line)) {
        ++line_no;
        if (line_no == 1) {
            header = line;
        }
        if (line_no == 2) {
            seq = line;
        }
        if (line_no == 4) {
            score = line;

            seq = mutate_seq (seq, score, scores_arr, phred_scores);
            out << header + "\n";
            out << seq + "\n";
            out << plus + "\n";
            out << score + "\n";
            line_no = 0;
        }
    }
    return 0;
}


