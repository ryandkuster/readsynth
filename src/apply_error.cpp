#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <vector>


char mutate_base (char base) {
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
    return rand() < prob * (RAND_MAX+1.0);
}


std::string mutate_seq (std::string seq, std::string score, double scores_arr[], std::string phred_scores) {
    std::string mut_seq = "";
    for (int i = 0; i < seq.length(); i++) {
        char base = seq[i];

        if (base == 'N') {
            mut_seq = mut_seq + base;
            continue;
        }

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
    std::string file_name = argv[1];
    std::string out_name = argv[2];

    std::ifstream file(file_name);
    std::ofstream out(out_name);

    std::string phred_scores = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ";
    int N = phred_scores.length();
    double scores_arr[N];

    for (int i = 0; i < N; i++) {
        double y = -(double)i/10;
        scores_arr[i] = pow(10, y);
    }

    int line_no = 0;
    std::string line;
    std::string header = "";
    std::string seq = "";
    std::string plus = "+";
    std::string score = "";

    srand(time(0));

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


