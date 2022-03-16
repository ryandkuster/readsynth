# readsynth

## what is readsynth?
Readsynth is a series of python scripts that simulate reduced-representation sequencing libraries with special consideration to DNA fragment abundance. Readsynth takes as input any reference nucleotide sequences (fasta format) along with custom restriction enzyme and adapter information to produce a simulatied set of reads expected in an Illumina-based sequencing reaction.

Readsynth is aimed at standardizing reduced-metagenome sequencing (RMS) libraries, where diverse community members are expected to yield differences in abundance due to molecular sequencing preparation.

## what makes readsynth different?
The preparation of DNA sequencing libraries include many variables which influence the frequency of a given read, including restriction enzyme site frequency in the sample genome(s), enzyme digestion efficiency, size-selection, as well as PCR bias based on DNA fragment length. Readsynth allows users to control common sources of variability introduced in typical enzymatic library preparations (e.g. RADseq, ddRADseq, GBS, etc.).

## requirements
Python packages numpy, pandas, and seaborn
To install using pip:
```
python3 -m pip install numpy
python3 -m pip install pandas
python3 -m pip install seaborn
```

To install using conda (in a conda environment):
```
conda install numpy
conda install pandas
conda install seaborn
```

Command line core utilities, specifically 'shuf'. To see if shuf is available on your system, try:
```
shuf --help
```

To install coreutils on macOS:
```
brew install coreutils
```

To install coreutils using the apt package manager:
```
sudo apt-get -y install coreutils
```

## usage

**inputs**
- genome assembly file (fasta) OR csv of locally stored genome files
- optional: custom adapter sequences
- optional: pre-existing fastq data to train error profile
- see full list of custom settings under 'input options' below

**outputs**
- csv of **all possible fragments** within a specified maximum length (base pairs)
- csv of **expected fragment counts** given a specified restriction enzyme digestion efficiency
- csv of **size-selected fragment counts** within a specified Normal distribution of read lengths
- optional: simulated fastq file of expected reads


## example
```
readsynth.py -genome reference.fasta -m1 G/AATTC -m2 T/TAA -n 10000 -cut_prob 0.90 -mean 300 -sd 100 -l 100 -o /output_directory
```
The above example takes 'reference.fasta' as a **genome** to be digested with EcoRI (**m1**) and MseI (**m2**) in a strand-specific fashion (e.g. the forward adapter ligates with the /AATTC overhang from EcoRI). Assuming 10,000 genomic copies (**n**) of the genome, digest simulation will calculate the expected number of DNA fragments given the enzyme digestion efficiency (**cut_prob**) occurs at a probability of 90% at any random RE motif. The resulting fragments will be size-selected using a normal distribution defined by **mean** and **sd**. Paired-end Illumina reads of length (**l**) will be written to a simulated fastq file (default output has perfect scores).

## input options
- genome  - path to file genome
- o - path to store output
- m1 - space separated list of RE motifs (e.g., AluI = AG/CT,
                    HindIII = A/AGCTT, SmlI = C/TYRAG)
- m2 - space separated list of RE motifs (e.g., AluI = AG/CT,
                    HindIII = A/AGCTT, SmlI = C/TYRAG)
- l - desired read length of final simulated reads (defaults
                    to 250 or given q1/q2 profiles)
- test - test mode: create newline-separated file of RE digested
                    sequences only
- t - number of subprocesses to run while simulating copy
                    number
- n - genome copy (depth per locus)
- mean - mean (in bp) of read lengths after size selection
- sd - standard deviation (in bp) of read lengths after size
                    selection
- min - min distance between cuts (optional, defaults to 6bp)
- max - max fragment length after first cut (optional, defaults
                    to mean + 6 stdevs)
- cut_prob - percent probability of per-site cut; use '1' for
                    complete digestion of fragments (fragments will not
                    contain internal RE sites)
- a1 - file containing tab/space-separated adapters and barcode
                    that attach 5' to read
- a2 - file containing tab/space-separated adapters and barcode
                    that attach 3' to read
- a1s - manually provide bp length of adapter a1 before SBS
                    begins
- a2s - manually provide bp length of adapter a1 before SBS
                    begins
- q1 - file containing R1 q scores in csv format (see
                    ngsComposer tool crinoid)
- q2 - file containing R2 q scores in csv format (see
                    ngsComposer tool crinoid)
- r1 - R1 fastq file to sample Q scores
- r2 - R2 fastq file to sample Q scores
- p - if using r1/r2 for profile, percent of reads to sample

# software overview

```
--readsynth.py
 |||
 || `_ 1. digest_genomes.py  -  store seq and pos of possible RE fragments
 | `__ 2. n_copies.py        -  partially digest genome with n copies
  `___ 3. size_selection.py  -  size select reads from Gaussian distribution
```
## 1. digest_genomes.py  
Given a fastq sequence, digest_genomes.py performs a regex search for all m1/m2 restriction motifs. Once a motif hit occurs, the sequence is searched forward within a set range (**max**, default is mean + 6sd). The start and end position as well as the orientation of the motifs is saved, and sequences are reverse-complemented when oriented on the reverse direction of the reference sequence. Only fragments that contain the m1 motif(s) at the 5' end and the m2 motif(s) at the 3' end will be kept.

Below is a toy genome example of a search for m1=G/AATTC and m2=T/TAA where the maximum fragment length to search forward is 50bp:
```
     _______maximum expected fragment length__________
    |                                                 |

    4                                37         48                     71
    m1                               m2         m2                     m1
    v                                v          v                      v
GTGAGAATTCGTTGAAAATCCGGTCCTGACGGGACTTTTAACAAGGAATTAAAGATCGCCATAATATTATTGAATTCCC
CACTCTTAAGCAACTTTTAGGCCAGGACTGCCCTGAAAATTGTTCCTTAATTTCTAGCGGTATTATAATAACTTAAGGG

possible fragments:

     AATTCGTTGAAAATCCGGTCCTGACGGGACTTT
         GCAACTTTTAGGCCAGGACTGCCCTGAAAAT

     AATTCGTTGAAAATCCGGTCCTGACGGGACTTTTAACAAGGAAT
         GCAACTTTTAGGCCAGGACTGCCCTGAAAATTGTTCCTTAAT

                                      TAACAAGGAAT      m2 -> m2 orientation, fragment not kept
                                        TGTTCCTTAAT

                                      TAACAAGGAATTAAAGATCGCCATAATATTATTG
                                        TGTTCCTTAATTTCTAGCGGTATTATAATAACTTAA

                                                 TAAAGATCGCCATAATATTATTG
                                                   TTCTAGCGGTATTATAATAACTTAA
```
The resulting fragments are saved in the output directory as a csv-formatted 'raw_digest' file with the following information:


|seq|start|end|m1|m2|length|reverse|
|:- |:-   |:- |:-|:-|:-    |:-     |
|AATTCGTTGAAAATCCGGTCCTGACGGGACTTT|4|37|GAATTC|TTAA|37|0|
|AATTCGTTGAAAATCCGGTCCTGACGGGACTTTTAACAAGGAAT|4|48|GAATTC|TTAA|48|0|
|AATTCAATAATATTATGGCGATCTTTAATTCCTTGT|37|71|TTAA|GAATTC|40|1|
|AATTCAATAATATTATGGCGATCTT|48|71|TTAA|GAATTC|29|1|

Above, 'length' reflects the overhang-to-overhang distance and the sequences ('seq') are represented with the overhangs cleaved in the 5' to 3' orientation.

## 2. n_copies.py
Once the raw digest fragments have been determined, n_copies.py uses the per-site probability of enzyme cleaving (**cut_prob**) to simulate digesting **n** copy numbers of the input genome. This model captures the non-uniform distribution of fragments where fragments containing intervening cut sites will be less abundant.

Below is a sample output of our toy example using cut_prob = 90% and n = 10,000 genome copies:

|seq|start|end|m1|m2|length|reverse|copies|
|:- |:-   |:- |:-|:-|:-    |:-     |:-    |
|AATTCGTTGAAAATCCGGTCCTGACGGGACTTT|4|37|GAATTC|TTAA|37|0|8062|
|AATTCGTTGAAAATCCGGTCCTGACGGGACTTTTAACAAGGAAT|4|48|GAATTC|TTAA|48|0|805|
|AATTCAATAATATTATGGCGATCTTTAATTCCTTGT|37|71|TTAA|GAATTC|40|1|794|
|AATTCAATAATATTATGGCGATCTT|48|71|TTAA|GAATTC|29|1|8132|

Notice the second and third sequences are less abundant as they are mutually exclusive within a single genome (e.g. if cut_prob = 100%, only those fragments with no internal enzyme motifs will be produced in simulation).

## 3. size_selection.py
Now that we have a baseline number of fragments expected from n number of genome copies with a given digestion efficiency, we use the **mean** and **sd** arguments to draw fragments from a Normal distribution. The area under the distribution (number of reads) is grown until reads in the range of mean through mean + 2sd are not under-represented in the sampling process.

## fastq output formatting
In each sequence header, the following components may be useful for a number of analyses:

```
@62:47:138:34:CGAAGGTGAT:TGCCAAT:full_genome.fna 1
```

|@|62|47|138|34|CGAAGGTGAT|TGCCAAT|full_genome.fna|1|
|:-|:-|:-|:-|:-|:---------|:------|:-------------|:-|
|'@'|fastq index|index of 'sampled_df.csv'|length of fragment with adapters|length of fragment|R1 barcode|R2 barcode|fasta name|r1/r2 read|
