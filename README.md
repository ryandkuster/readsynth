readsynth is currently in an early release state. If you have any suggestions, questions, or concerns, please open an issue or contact rkuster@utk.edu.

<img src=resources/images/readsynth_logo_blue.png width="600">

## what is readsynth?
Readsynth simulates Illumina sequenced metagenomes allowing comparison of library approaches including amplicon and reduced metagenome sequencing (RMS). The software takes as input any reference genome sequences in fasta format along with custom primer/restriction enzyme and adapter sequences to simulate the reads expected from an Illumina-based sequencing reaction. Simulated fastq reads capture the expected composition of the input taxa, and are sensitive to changes in the underlying distributions of DNA fragments when digests or size selection are performed. 

Readsynth is aimed at optimizing sequencing effort and library prep settings to accurately profile diverse communities of microbial taxa. Use cases may include:
- determining sequencing effort (total reads) required to capture rare taxa
- selecting primers or restriction enzyme schema to capture species richness and diversity
- creating a ground truth for benchmarking profiling software

## contents
- [overview](#overview)
- [usage](#usage)
  - [installation](#installation)
  - [requirements](#requirements)
  - [memory considerations](#memory-considerations)
  - [example ddRADseq library simulation](#example-ddRADseq-library-simulation)
  - [example amplicon library simulation](#example-amplicon-library-simulation)
  - [example isolength enzyme library simulation](#example-isolength-enzyme-library-simulation)
  - [optional parameters](#optional-parameters)
- [software details](#software-details)
- [license](#license)

## overview
The preparation of DNA sequencing libraries includes many factors which influence the frequency of a given read, including restriction enzyme site frequency (taxa-dependent), enzyme digestion efficiency, size-selection, as well as the resulting biases caused by preferential DNA fragment length. Readsynth allows users to control common sources of variability introduced in typical enzymatic library preparations.

**inputs**
- table of locally stored genome assembly files (fasta) and abundance
- desired total number of reads
- restriction enzyme motifs and cut efficiency (can use iso-length REs or 16S/ITS primers)
- size distribution parameters (or locally stored json of custom distribution)
- optional: custom adapter sequences
- optional: pre-existing fastq data to train error profile
- see full list of custom settings in 'optional parameters'

<img src=resources/images/abundances.png width="600">

**outputs**
- simulated fastq file of expected reads (paired end)
- csv of **all possible fragments** within a specified maximum length (base pairs)
- csv of **expected fragment counts** given a specified restriction enzyme digestion efficiency
- csv of **size-selected fragment counts** within a specified Normal distribution of read lengths
- csv summary of reads produced
- individual and stacked distributions before and after size selection (see images below)

<img src=resources/images/_fragment_distributions.png width="600">
<img src=resources/images/_read_distributions.png width="600">

## usage

### installation

First, download readsynth from this github page.

To begin using readsynth for the first time, change to "src" within the readsynth directory and run:
```
make apply_error
```

### requirements
Python 3.7 and higher

Python packages numpy, pandas, and seaborn

To install using pip:
```
python3 -m pip install numpy
python3 -m pip install pandas
python3 -m pip install seaborn
```

Optionally, to install using conda (within a single conda environment):
```
conda install numpy
conda install pandas
conda install seaborn
```

### memory considerations

Readsynth is designed for assessing bacterial and fungal genomes and currently isn't optimized for use with larger eukaryotic genomes. To include larger genomes, it is recommended to split these fasta files into individual sequence files and input each of those into the abundance table. Memory usage will vary: a large sequence length combined with highly abundant motif sites will require enough ram to store every possible fragment produced in the expected fragment distribution.

### example ddRADseq library simulation

```
python3 readsynth.py -g abundances.csv -m1 EcoRI -m2 T/TAA -n 1_000_000 -c 0.90 -u 300 -sd 100 -l 150 -o /output_directory
```

The above example takes 'abundances.csv' as a genome abundance file **g** with all reference fasta files to be digested with EcoRI (**m1**) and MseI (**m2**) in a strand-specific fashion (e.g. the forward adapter always ligates with the /AATTC overhang from EcoRI). Assuming a desired output of 1 million reads (**n**) will be approximated, digest simulation will calculate the expected number of DNA fragments given the enzyme digestion cut efficiency (**c**) occurs at a probability of 90% at any random RE motif. The resulting fragments will be size-selected using a normal distribution defined by **u** and **sd**. Paired-end Illumina reads of length (**l**) 150bp will be written to a simulated fastq file (default output has perfect scores). All output (**o**) will be written to './output_directory' relative to the present working directory where readsynth is called.

<sub>relevant input options for reduced sequencing library:</sub>
```
  -g G             path to file genome
  -o O             path to store output
  -m1 M1 [M1 ...]  space separated list of search motifs (e.g., RE motifs AluI or AG/CT)
  -m2 M2 [M2 ...]  space separated list of search motifs (e.g., RE motifs SmlI or C/TYRAG)
  -l L, -l1 L      desired R1 read length of final simulated reads
  -n N             total read number
  -u U             mean (in bp) of read lengths after size selection
  -sd SD           standard deviation (in bp) of read lengths after size selection
  -x X             fragment length where fragment distribution intersects size distribution
  -d D             json dictionary of fragment length:count for all expected bp fragments range
  -lp LP           low-pass mode: defines maximum expected fragment size, distribution free
```

### example amplicon library simulation

```
python3 readsynth.py -g abundances.csv -m1 /CCTACGGGNGGCWGCAG -m2 /GACTACHVGGGTATCTAANCC -n 1_000_000 -lp 1000 -l 150 -o /output_directory
```

The above example takes 'abundances.csv' as a genome abundance file **g** with all reference fasta files to be PCR targeted with forward (**m1**) and reverse (**m2**) 16S primer sequences in a strand-specific fashion (e.g. the forward adapter always ligates with the forward primer). A desired output of 1 million reads (**n**) will be approximated. Further, this library avoids applying a gaussian size selection step and utilizes the **lp** argument to include any amplicon sequences up to a maximum "low-pass" length of 1000bp. Paired-end Illumina reads of length (**l**) 150bp will be written to a simulated fastq file (default output has perfect scores). All output (**o**) will be written to './output_directory' relative to the present working directory where readsynth is called.

<sub>relevant input options for amplicon library:</sub>
```
  -g G             path to file genome
  -o O             path to store output
  -m1 M1 [M1 ...]  space separated list of forward primers (e.g., /CCTACGGGNGGCWGCAG)
  -m2 M2 [M2 ...]  space separated list of forward primers (e.g., /GACTACHVGGGTATCTAANCC)
  -l L, -l1 L      desired R1 read length of final simulated reads
  -n N             total read number
  -lp LP           low-pass mode: defines maximum expected amplicon size
```

### example isolength enzyme library simulation

```
python3 readsynth.py -g abundances.csv -iso BcgI -n 1_000_000 -l 150 -lp 50 -o /output_directory
```

The above example takes 'abundances.csv' as a genome abundance file **g** with all reference fasta files to be digested with the isolength (**iso**) type IIb enzyme BcgI. A desired output of 1 million reads (**n**) will be approximated. Further, this library avoids applying a gaussian size selection step and utilizes the **lp** argument to include fragments up to a maximum "low-pass" length of 50bp. Paired-end Illumina reads of length (**l**) 150bp will be written to a simulated fastq file (default output has perfect scores). All output (**o**) will be written to './output_directory' relative to the present working directory where readsynth is called.

<sub>relevant input options for isolength library:</sub>
```
  -g G             path to file genome
  -o O             path to store output
  -iso ISO         optional type IIB RE motif (e.g., BcgI or NN/NNNNNNNNNNCGANNNNNNTGCNNNNNNNNNNNN/)
  -l L, -l1 L      desired R1 read length of final simulated reads
  -n N             total read number
  -lp LP           low-pass mode: defines maximum expected fragment size
```

### optional parameters

```
  -c C             optional: percent probability of per-site cut; default 1 for complete digestion of fragments (fragments will not contain internal RE sites)
  -a1 A1           optional: file containing tab/space-separated adapters and barcode that attach 5' to read
  -a2 A2           optional: file containing tab/space-separated adapters and barcode that attach 3' to read
  -a1s A1S         optional: manually provide bp length of adapter a1 before SBS begins
  -a2s A2S         optional: manually provide bp length of adapter a1 before SBS begins
  -q1 Q1           optional: file containing newline-separated R1 Q scores >= length -l
  -q2 Q2           optional: file containing newline-separated R2 Q scores >= length -l
  -e E             optional: filler base to use if full adapter contaminaton occurs
  -l2 L2           optional: desired R2 read length of final simulated reads
  -test            test mode: skip writing simulated fastq files
```


## software details

```
--readsynth.py
|| |
|| |`_ 1a. digest_genomes.py -  store seq and pos of possible RE fragments
|| `__ 1b. prob_n_copies.py  -  partially digest genome with n copies
| `___ 2.  size_selection.py  -  size select reads from defined distribution
 `____ 3.  write_reads.py     -  write fragments as paired end fastq-formatted reads
```

### 1. digest_genomes.py  
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
The resulting fragments are saved in the output directory as a csv-formatted 'raw_digest' file.

### 2. prob_n_copies.py
Once the raw digest fragments have been determined, prob_n_copies.py uses the per-site probability of enzyme cleaving (**c**) to simulate digesting **n** copy numbers of the input genome. This model captures the non-uniform distribution of fragments where fragments containing intervening cut sites will be less abundant.

### 3. size_selection.py
Now that we have a baseline number of fragments expected from n number of genome copies with a given digestion efficiency, we use the **mean** and **sd** arguments to draw fragments from a Normal distribution. The area under the distribution (number of reads) is grown until reads in the range of mean through mean + 2sd are not under-represented in the sampling process.

### fastq output formatting
In each sequence header, the following components may be useful for a number of analyses:

```
@62:47:138:34:CGAAGGTGAT:TGCCAAT:full_genome.fna 1
```

|@|62|47|138|34|CGAAGGTGAT|TGCCAAT|full_genome.fna|1|
|:-|:-|:-|:-|:-|:---------|:------|:-------------|:-|
|'@'|fastq index|index of counts file|length of fragment with adapters|length of fragment|p5 adapter id|p7 adapter id|fasta name|r1/r2 read|

## license

<a href="https://github.com/ryandkuster/readsynth/blob/master/LICENSE">Apache-2.0 license</a>
