---
title: "Finding and downloading raw data"
teaching: 0
exercises: 0
questions:
- "Key question (FIXME)"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

In this tutorial we will be looking at some published data on
the binding of Gata1 and Gata2 in G1E cells (an erythroid cell line
in which Gata1 is fused to an ER domain, allowing for induction of
terminal erythroid differentiation using tamoxifen).

The paper with this information is [Trompouki et al.][g1e-paper],
published in Cell in 2011.

SRX100313

We need to download the file from the EBI direct to our server. We can do this
with a programme called wget:

~~~
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR351/SRR351407/SRR351407.fastq.gz
~~~
{: .language-bash}

This has downloaded one file called SRR351407.fastq.gz

Let's try to examine the contents of the file:

~~~
head SRR351407.fastq.gz
~~~
{: .language-bash}

~~~
��}��85��"U�<�3c��N@�H����ĩJX�|�O�O������3����	�~9��������m����������?�
o��������3թ>��R�I��{�������|�u]�V�ժmZoUg���MW{��qVjߵǮp+Xmt�y��qV;�j_��[U�Y
�%��9Ϊ���j�}V-Z�ű{�,$����}�R�N�k���J�Iw�`�TJ��,��:v#���*˲Uە|�7�Ux��@Ò�*˾Q�.
z�RrTe�V��,�R}y�a�]�u��*�_��R�X�Q��TS6�X�QI�l�����]V�TG�]VӐ�:r_�.�i���T�.�
鋣}e%����+�eu�.~â�����U������<�=V�����+z�^�G?��T���
~~~
{: .output}

We get a lot of nonsense because this is a zipped file. We can use a command
called zcat to unzip these types of files and print them to the terminal.
Since this file has millions of lines, we want to make sure we pipe the output
to head.

~~~
zcat SRR351407.fastq.gz | head -n 12
~~~
{: .language-bash}

~~~
@SRR351407.1 WICMT-SOLEXA2:3:1:1801:997/1
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
+
####################################
@SRR351407.2 WICMT-SOLEXA2:3:1:2545:999/1
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
+
####################################
@SRR351407.3 WICMT-SOLEXA2:3:1:3463:995/1
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
+
####################################
~~~
{: .output}

A FASTQ file normally uses four lines per sequence.

- Line 1 begins with a '@' character and is followed by the read name
- Line 2 is the raw DNA sequence
- Line 3 begins with a '+' character
- Line 4 encodes the quality values for the sequence in Line 2, and must contain the same number of symbols as letters in the sequence.

So our first read is called `SRR351407.1 WICMT-SOLEXA2:3:1:1801:997/1` and its sequence
is `NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN`. We obtained no useful sequence information from
this read. This is quite common at the beginning of a fastq file, especially from
older Illumina machines. We need to look at some reads which are further into the file.

> ## Peeking inside the middle of a file
>
> Write a command that will allow us to view lines
> 10000-10012 of SRR351407.fastq.gz
>
> > ## Solution
> >
> > You need to combine zcat, tail and head in a pipe to do this
> >
> > ~~~
> > zcat SRR351407.fastq.gz | head -n 10012 | tail -n 12
> > ~~~
> > {: .language-bash}
> >
> > ~~~
> > @SRR351407.2501 WICMT-SOLEXA2:3:1:7736:1233/1
> > ACACCTTTTCCTGCAGGGACATCGTCTGCCACCGAC
> > +
> > GGGEGGEDEDGGGGG>3?5AG@GGDGGEGD8BAAAD
> > @SRR351407.2502 WICMT-SOLEXA2:3:1:7847:1231/1
> > GTCACTGCCACTGTTGGCCAGGATGCACACACACAC
> > +
> > IIIIIIIIIIIIIIGIIIIFIIEIIEEGEG@IGGGG
> > @SRR351407.2503 WICMT-SOLEXA2:3:1:7925:1234/1
> > GGAGGCAAACCTGACCTGCCTTCCCTGTAACGGTGG
> > +
> > GIIIIIIHIIIIIIIIIIIIIIIIIIIHIIFIIDID
> > ~~~
> > {: .language-bash}
> {: .solution}
{: .challenge}

These fastq files are big, and will take a long time to process. To
be able to complete this tutorial in a reasonable amount of time, we
can download some data files which have had most of the reads removed,
leaving only those that map to chromosome 11. We use `wget` for this,
and then unzip the files using a program called `tar`.

~~~
wget https://rob.beagrie.com/media/chip-tutorial/chip-tutorial-files.tar.gz
tar zxvf chip-tutorial-files.tar.gz
cd chip-tutorial
ls fastqs
~~~
{: .language-bash}

~~~
G1E_ATAC_0h_r1A_end1.chr11.fastq.gz  G1E_ATAC_30h_r1_end1.chr11.fastq.gz
G1E_ATAC_0h_r1A_end2.chr11.fastq.gz  G1E_ATAC_30h_r1_end2.chr11.fastq.gz
G1E_ATAC_0h_r1B_end1.chr11.fastq.gz  G1E_ATAC_30h_r2_end1.chr11.fastq.gz
G1E_ATAC_0h_r1B_end2.chr11.fastq.gz  G1E_ATAC_30h_r2_end2.chr11.fastq.gz
G1E_ATAC_0h_r1_end1.chr11.fastq.gz   SRR351406.chr11.fastq.gz
G1E_ATAC_0h_r1_end2.chr11.fastq.gz   SRR351407.chr11.fastq.gz
G1E_ATAC_0h_r2_end1.chr11.fastq.gz   SRR351409.chr11.fastq.gz
G1E_ATAC_0h_r2_end2.chr11.fastq.gz   SRR351410.chr11.fastq.gz
~~~
{: .output}

Our example data contains some ATAC-seq data files, which have
sensible names. We've also downloaded four files from
[Trompouki et al.][g1e-paper] which are named using their ID
from GEO. It's going to be difficlut to remember what these cells are.
We could rename them, but then we might forget which file relates to
which piece of published data. One solution here is to make a symbolic
link using `ln -s`

~~~
cd fastqs
ln -s SRR351406.chr11.fastq.gz G1E_ChIP_Input_0h.chr11.fastq.gz
ln -s SRR351407.chr11.fastq.gz G1E_ChIP_Gata2_0h.chr11.fastq.gz
zcat G1E_ChIP_Input_0h.chr11.fastq.gz | head
~~~
{: .language-bash}

> ## Linking files
>
> Make symbolic links to name the two remaining files from the
> G1E paper following the same convention.
>
> > ## Solution
> >
> > SRR351409 and SRR351410 are the input and the Gata1 ChIP
> > from G1ER cells treated for 24h with estradiol. So a
> > sensible name would look something like this:
> >
> > ~~~
> > ln -s SRR351409.chr11.fastq.gz G1E_ChIP_Input_24h.chr11.fastq.gz
> > ln -s SRR351410.chr11.fastq.gz G1E_ChIP_Gata1_24h.chr11.fastq.gz
> > ~~~
> > {: .language-bash}
> {: .solution}
{: .challenge}



{% include links.md %}

