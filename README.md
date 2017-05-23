# HiRGC

A high performance referential genome compression algorithm (termed HiRGC)

#### Result
The performance of HiRGC is compared with four state-of-the-art compression algorithms on eight human genome benchmark data sets. HiRGC uses less than 30 minutes to compress about 21 gigabytes of each set of seven target genomes into 96 to 260 megabytes, achieving compression ratios of 217 to 82 times. This performance is one order of magnitude better than most of the state-of-the-art algorithms, and it is at least 1.9 times better than their best performance. The compression speed is at least 2.9 times faster. HiRGC also exhibits a very stable and robust performance when tested on different reference genomes, greatly narrowing down the wide variation of the existing algorithms’ performance.

Two genomes (i.e., HG19 and YH) and a compressed result by iDoComp can be downloaded from [here](https://www.dropbox.com/s/3lg2131t2wdo6la/hg19_YH.zip?dl=0 "genomes HG19 and YH")

## Compile
	make hirgc
	make de_hirgc

## Compress
The following three different commands:

    （1） ./hirgc -m file -r YH_chr1.fa -t HG18_chr1.fa
    （2） ./hirgc -m genome -r YH -t HG18 -n chr_name.txt/default
    （3） ./hirgc -m set -r YH -t genome_set.txt -n chr_name.txt/default


* -m is the mode that the user want to used;
* -r is the reference (a FASTA file or a genome folder);
* -t is the target (a FASTA file or a genome folder or a file contains a list of genome folder);
* -n is a file containing name of chromosomes or a string "default".

---
**Some other explaination**<br />

- the mode *file* is to comrpess chromosomes; *genome* is to comoress a genome; *set* is to compress genome set.
- The first command compress a chromosome *HG18_chr1.fa* and the reference chromosome is *YH_chr1.fa*. The compressed is *HG18_chr1.fa_ref_YH_chr1.fa.7z*.
- The second compress a genome *HG18* contain some chromosomes and the reference genome is *YH*. The name of chromosomes can be given in the file *chr_name.txt* or set as default. The compressed file is *HG18_ref_YH.7z*.
- The third compress a genome set and the reference genome is *YH*. The list of genomes is given in file *genome_set.txt*. The name of chromosomes can be given in the file *chr_name.txt* or can be set as default.
- The default name of chromosomes is [*chr1.fa, chr2.fa, chr3.fa, chr4.fa, chr5.fa, chr6.fa, chr7.fa, chr8.fa, chr9.fa, chr10.fa, chr11.fa, chr12.fa, chr13.fa, chr14.fa, chr15.fa, chr16.fa, chr17.fa, chr18.fa, chr19.fa, chr20.fa, chr21.fa, chr22.fa, chrX.fa, chrY.fa*]

## Decompress

    ./de_hirgc -m file -r YH_chr1.fa -t HG18_chr1.fa_ref_YH_chr1.fa.7z
    ./de_hirgc -m genome -r YH -t HG18_ref_YH.7z -n chr_name.txt/default
    ./de_hirgc -m set -r YH -t de_genome_set.txt -n chr_name.txt/default

- the first command decompress a compressed chromosome file *xxx.7z* to *dec_xxx*.
- the second command decomperss a  compressed genome folder *yyy.7z* to *dec_yyy*.
- the parameter *-n* is identical to the compression procedure.

## Status
Submitted to *Bioinformatics* (major revision).

### Contacts
If any bugs during you run our code, please email to <yyuanshengliu@gmail.com>

