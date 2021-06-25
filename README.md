## readsynth

### code architecture
what should the system do? (functional requirements)
- digest input genomes or chromosomes into RAD libraries

how should the system behave? (non-functional requirements)

restrictions
- timeframe
- existing tools


### simulated reads formatting

in each header, the following components may be useful for analysis:

```
@62:47:138:34:CGAAGGTGAT:TGCCAAT:full_genome.fna 1
```

|@|62|47|138|34|CGAAGGTGAT|TGCCAAT|full_genome.fna|1|
|:-|:-|:-|:-|:-|:---------|:------|:-------------|:-|
|'@'|fastq index|index of 'sampled_df.csv'|length of fragment with adapters|length of fragment|R1 barcode|R2 barcode|fasta name|r1/r2 read|

