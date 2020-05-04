---
title: "Aligning data to a genome"
teaching: 0
exercises: 0
questions:
- "Key question (FIXME)"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

We now need to figure out where these reads came from in the
mouse genome. There are a class of programs specially designed
for this task, generally called "aligners". The one we will be
using for this tutorial is `bowtie2`. On the cluster, we need
to specify which version of bowtie2 we want to use before
we can start aligning our reads. To see the list of available
versions, we can use `module avail`

Once bowtie2 is loaded, we can then view its help page to give
us an idea what to do next.

~~~
module avail
module load bowtie2/2.3.5
bowtie2 --help
~~~
{: .language-bash}

~~~
Bowtie 2 version 2.3.5 by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea)
Usage:
  bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | --sra-acc <acc>} [-S <sam>]

  <bt2-idx>  Index filename prefix (minus trailing .X.bt2).
             NOTE: Bowtie 1 and Bowtie 2 indexes are not compatible.
  <m1>       Files with #1 mates, paired with files in <m2>.
             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
  <m2>       Files with #2 mates, paired with files in <m1>.
             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
  <r>        Files with unpaired reads.
             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
  <i>        Files with interleaved paired-end FASTQ/FASTA reads
             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
  <acc>      Files are SRA accessions. Accessions not found in local storage will
             be fetched from NCBI.
  <sam>      File for SAM output (default: stdout)

  <m1>, <m2>, <r> can be comma-separated lists (no whitespace) and can be
  specified many times.  E.g. '-U file1.fq,file2.fq -U file3.fq'.
~~~
{: .output}

The first thing we need to provide bowtie is an index. Normally, we would use
an index that covers the whole mouse genome. For this tutorial, we can use an
index that only covers mouse chromosome 11. This will make bowtie2 run more
quickly and use less memory.

~~~
ls bowtie2-index-mm10-chr11
~~~
{: .language-bash}

~~~
genome.1.bt2  genome.2.bt2  genome.3.bt2
genome.4.bt2  genome.rev.1.bt2  genome.rev.2.bt2
~~~
{: .output}

The other thing we need to provide to bowtie is a fastq file to map.
The files downloaded from [Trompouki et al.][g1e-paper] are "unpaired"
or "single-end" reads, so we need to use the `-u` flag.

~~~
bowtie2 -x bowtie2-index-mm10-chr11/genome -U fastqs/G1E_ChIP_Gata1_24h.chr11.fastq.gz | head -n 4
~~~
{: .language-bash}

~~~
@HD	VN:1.0	SO:unsorted
@SQ	SN:chr11	LN:122082543
@PG	ID:bowtie2	PN:bowtie2	VN:2.3.5	CL:"/package/bowtie2/2.3.5/bin/bowtie2-align-s --wrapper basic-0 -x bowtie2-index-mm10-chr11/genome -U fastqs/G1E_ChIP_Gata1_24h.chr11.fastq.gz"
SRR351410.24963042	16	chr11	3102483	0	36M	*	0	0	GCATTTCACATTATTCACGTTTTTCACTGTTTCTCG	FFFFFFFFFFFGGEGGGGGGFGGGGGGFGGGGFGFG	AS:i:-10	XS:i:-10	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:1T10T23	YT:Z:UU
~~~
{: .output}

The first three lines start with @ - these are the header lines that give information
about the command that was used to generate the file and the index that was used
for the mapping.

The last line is a single alignment of a read agains the genome. This information is given
in SAM format:

| Col     | Field Type | Brief description                     | First read of our   file             |
|---------|------------|---------------------------------------|--------------------------------------|
| 1 QNAME | String     | Query template NAME                   | SRR351410.24963042                   |
| 2 FLAG  | Int        | bitwise FLAG                          | 16                                   |
| 3 RNAME | String     | Reference sequence NAME               | chr11                                |
| 4 POS   | Int        | 1-based leftmost mapping POSition     | 3102483                              |
| 5 MAPQ  | Int        | MAPping Quality                       | 0                                    |
| 6 CIGAR | String     | CIGAR string                          | 36M                                  |
| 7 RNEXT | String     | Reference name of the mate/next read  | *                                    |
| 8 PNEXT | Int        | Position of the mate/next read        | 0                                    |
| 9 TLEN  | Int        | observed Template LENgth              | 0                                    |
| 10 SEQ  | String     | segment SEQuence                      | GCATTTCACATTATTCACGTTTTTCACTGTTTCTCG |
| 11 QUAL | String     | ASCII of Phred-scaled base QUALity+33 | FFFFFFFFFFFGGEGGGGGGFGGGGGGFGGGGFGFG |

Let's have a look at another of our files:

~~~
bowtie2 -x bowtie2-index-mm10-chr11/genome -U fastqs/G1E_ChIP_Gata2_0h.chr11.fastq.gz | head -n 4
~~~
{: .language-bash}

~~~
@HD	VN:1.0	SO:unsorted
@SQ	SN:chr11	LN:122082543
@PG	ID:bowtie2	PN:bowtie2	VN:2.3.5	CL:"/package/bowtie2/2.3.5/bin/bowtie2-align-s --wrapper basic-0 -x bowtie2-index-mm10-chr11/genome -U fastqs/G1E_ChIP_Gata2_0h.chr11.fastq.gz"
SRR351406.214265	16	chr11	3099998	0	12M1I23M	*	0	0	TCCAGGTCTTACGGTGTGTATTTCACATTTTTCACG	DGGDGG<GGEDDGGGGGBG8GGGG<GGGGGGGGGGG	AS:i:-21	XN:i:3	XM:i:5	XO:i:1	XG:i:1	NM:i:6	MD:Z:0N0N0N5C3T22	YT:Z:UU
~~~
{: .output}

Note the different CIGAR string: "12M1I23M" indicating that our read had one base more than
the reference genome at position 13.

We're going to need to save these alignments for future use:

~~~
mkdir aligned_reads
bowtie2 -x bowtie2-index-mm10-chr11/genome -U fastqs/G1E_ChIP_Gata2_0h.chr11.fastq.gz > aligned_reads/G1E_ChIP_Gata2_0h.chr11.sam
ls -l -h fastqs/SRR351407.chr11.fastq.gz aligned_reads/G1E_ChIP_Gata2_0h.chr11.sam
~~~
{: .language-bash}

~~~
-rw-r--r-- 1 user group 98M Apr 30 17:21 aligned_reads/G1E_ChIP_Gata2_0h.chr11.sam
-rw-r--r-- 1 user group  15M Apr 28 10:43 fastqs/SRR351407.chr11.fastq.gz
~~~
{: .output}

We've generated a sam file with all the alignments, but it's huge! It's around
6x larger than the fastq file. This is because the fastq is compressed but the sam
file isn't. Generally, people store alignments in bam files rather than sam
files. These are essentially just a compressed version of a sam. We can create a
bam file from a sam file using `samtools`.

~~~
module load samtools/1.9
samtools view -b aligned_reads/G1E_ChIP_Gata2_0h.chr11.sam > aligned_reads/G1E_ChIP_Gata2_0h.chr11.bam
ls -l -h aligned_reads/G1E_ChIP_Gata2_0h.chr11.sam aligned_reads/G1E_ChIP_Gata2_0h.chr11.bam
~~~
{: .language-bash}

~~~
-rw-r--r-- 1 user group  20M May  1 09:42 aligned_reads/G1E_ChIP_Gata2_0h.chr11.bam
-rw-r--r-- 1 user group 98M May  1 09:41 aligned_reads/G1E_ChIP_Gata2_0h.chr11.sam
~~~
{: .output}

That's better, the bam file is only slightly larger than the fastq. We
can do this in a single command so that we don't need to generate the
intermediate sam file.

~~~
rm aligned_reads/*
bowtie2 -x bowtie2-index-mm10-chr11/genome -U fastqs/G1E_ChIP_Gata2_0h.chr11.fastq.gz | samtools view -b > aligned_reads/G1E_ChIP_Gata2_0h.chr11.bam
ls -l -h aligned_reads
~~~
{: .language-bash}

~~~
-rw-r--r-- 1 user group  20M May  1 09:42 aligned_reads/G1E_ChIP_Gata2_0h.chr11.bam
~~~
{: .output}


> ## Mapping files in a for loop
>
> Using the command we just ran, write a for loop that will
> align all the files matching the pattern fastqs/G1E_ChIP_*.fastq.gz
>
> Hint: You will need to generate the bam files in fastqs/ and then
> move them to aligned_reads/
>
> > ## Solution
> >
> > ~~~
> > for fastq in fastqs/G1E_ChIP_*.fastq.gz;
> > do
> >    bowtie2 -x bowtie2-index-mm10-chr11/genome -U $fastq | samtools view -b > $fastq.bam
> >    mv $fastq.bam aligned_reads/
> > done
> > ~~~
> > {: .language-bash}
> >
> > or, without the move command:
> >
> > ~~~
> > cd fastqs
> > for fastq in G1E_ChIP_*.fastq.gz;
> > do
> >    bowtie2 -x ../bowtie2-index-mm10-chr11/genome -U $fastq | samtools view -b > ../aligned_reads/$(echo $fastq | cut -f 1,2 -d '.').bam
> > done
> > ~~~
> > {: .language-bash}
> {: .solution}
{: .challenge}

What about the ATAC-seq files? These are paired end data, so
there are two files for each dataset.

~~~
zcat fastqs/G1E_ATAC_0h_rep1_end1.chr11.fastq.gz | head -n 4
~~~
{: .language-bash}

~~~
@NB501183:592:HN75MBGX9:1:11101:1153:3424/1
TCTTCACTCTCCCTCTCTTCCGGATTCTCTTTAGGCCTTG
+
AAAAAEAEEEEEEEAEEAE/EEEEEEEEEEEEEE/EAEEE
~~~
{: .output}

~~~
zcat fastqs/G1E_ATAC_0h_rep1_end2.chr11.fastq.gz | head -n 4
~~~
{: .language-bash}

~~~
@NB501183:592:HN75MBGX9:1:11101:1153:3424/2
CACAAGGCCGAGAGCCCAGCCCGCAAGAAAGGCCCACGC
+
AAAAAAEEEEEEEE<E<<EE/EEEEEEEEEEEE<<AEEE
~~~
{: .output}

Note how the read name is the same, except that in the
`end1` file the read name ends in /1 whereas in the
`end2` file it ends in /2. These are the two "mates" or
"pairs" for this read.

To map these reads using bowtie we need to use the `-1` and
`-2` flags.

~~~
bowtie2 -x bowtie2-index-mm10-chr11/genome -1 fastqs/G1E_ATAC_0h_rep1_end1.chr11.fastq.gz -2 fastqs/G1E_ATAC_0h_rep1_end2.chr11.fastq.gz | samtools view -b > aligned_reads/G1E_ATAC_0h_rep1.bam
~~~
{: .language-bash}

What if we have more than one file for each experiment?
This is a common situation, especially with Illumina NextSeq machines
which give you four files for each sample. `bowtie2` needs us to give
these as a comma-separated list.

~~~
bowtie2 -x bowtie2-index-mm10-chr11/genome -1 fastqs/G1E_ATAC_0h_rep2A_end1.chr11.fastq.gz,fastqs/G1E_ATAC_0h_rep2B_end1.chr11.fastq.gz -2 fastqs/G1E_ATAC_0h_rep2A_end2.chr11.fastq.gz,fastqs/G1E_ATAC_0h_rep2B_end2.chr11.fastq.gz | samtools view -b > aligned_reads/G1E_ATAC_0h_rep2.bam
~~~
{: .language-bash}

> ## Finish mapping the files
>
> Make bam files of aligned reads for the rest of the atac-seq fastqs
>
> > ## Solution
> >
> > ~~~
> > bowtie2 -x bowtie2-index-mm10-chr11/genome -1 fastqs/G1E_ATAC_30h_rep1_end1.chr11.fastq.gz -2 fastqs/G1E_ATAC_30h_rep1_end2.chr11.fastq.gz | samtools view -b > aligned_reads/G1E_ATAC_30h_rep1.bam
> > bowtie2 -x bowtie2-index-mm10-chr11/genome -1 fastqs/G1E_ATAC_30h_rep2_end1.chr11.fastq.gz -2 fastqs/G1E_ATAC_30h_rep2_end2.chr11.fastq.gz | samtools view -b > aligned_reads/G1E_ATAC_30h_rep2.bam
> > ~~~
> > {: .language-bash}
> {: .solution}
{: .challenge}

> ## Mapping paired end files in a loop
>
> If we wanted to map the paired end files in a loop, the first thing to do would be
> to write a for loop that printed the correct reads and the correct bam file.
> Fill in the blanks to make the right loop:
>
> ~~~
> for timepoint in ___
> do
>   for replicate in rep1 rep2
>   ___
>     echo READS=fastqs/G1E_ATAC_${timepoint}____ fastqs/G1E_ATAC_${timepoint}____
>     echo BAM=___
>   done
> done
> ~~~
> {: .language-bash}
>
> ### Expected output:
>
> ~~~
> READS=fastqs/G1E_ATAC_0h_rep1_end1.chr11.fastq.gz fastqs/G1E_ATAC_0h_rep1_end2.chr11.fastq.gz
> BAM=aligned_reads/G1E_ATAC_0h_rep1.chr11.bam
> READS=fastqs/G1E_ATAC_0h_rep2_end1.chr11.fastq.gz fastqs/G1E_ATAC_0h_rep2_end2.chr11.fastq.gz
> BAM=aligned_reads/G1E_ATAC_0h_rep2.chr11.bam
> READS=fastqs/G1E_ATAC_30h_rep1_end1.chr11.fastq.gz fastqs/G1E_ATAC_30h_rep1_end2.chr11.fastq.gz
> BAM=aligned_reads/G1E_ATAC_30h_rep1.chr11.bam
> READS=fastqs/G1E_ATAC_30h_rep2_end1.chr11.fastq.gz fastqs/G1E_ATAC_30h_rep2_end2.chr11.fastq.gz
> BAM=aligned_reads/G1E_ATAC_30h_rep2.chr11.bam
> ~~~
> {: .output}
> > ## Solution
> >
> > ~~~
> > for timepoint in 0h 30h
> > do
> >   for replicate in rep1 rep2
> >   do
> >     echo READS=fastqs/G1E_ATAC_${timepoint}_${replicate}_end1.chr11.fastq.gz fastqs/G1E_ATAC_${timepoint}_${replicate}_end2.chr11.fastq.gz
> >     echo BAM=aligned_reads/G1E_ATAC_${timepoint}_${replicate}.chr11.bam
> >   done
> > done
> > ~~~
> > {: .language-bash}
> {: .solution}
{: .challenge}


{% include links.md %}

