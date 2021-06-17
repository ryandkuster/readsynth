## readsynth

### simulated reads formatting

in each header, the following components may be useful for analysis:

```
@9982:8912:CM008331.1_Ipomoea_batatas_cultivar_Taizhong6_chromosome_1,_whole_genome_shotgun_sequence:CTATGGATAA 1
```

|@|9982|8912|CTATGGATAA|AACCACTGG|CM008331.1_Ipomoea_batatas_...|1|
|:-|:--|:---|:---------|:---------|:-----------------------------|:-|
|'@'|final fastq index|sampled read index (matches 'sampled_df.csv')|R1 barcode|R2 barcode|fasta entry name|r1/r2 read|

