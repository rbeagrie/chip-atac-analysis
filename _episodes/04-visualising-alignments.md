---
title: "Visualising aligned reads"
teaching: 0
exercises: 0
questions:
- "Key question (FIXME)"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

To make sure we're all working with the same set of files, lets download
some bam files and overwrite the ones we made in the last lesson:

~~~
cd ..
rm -rf chip-tutorial/aligned_reads
wget https://rob.beagrie.com/media/chip-tutorial/chip-tutorial-files2.tar.gz
tar zxvf chip-tutorial-files2.tar.gz
ls -l chip-tutorial/aligned_reads
~~~
{: .language-bash}

~~~
-rw-r--r-- 1 user group  84893348 May  5 15:27 G1E_ATAC_0h_rep1.chr11.bam
-rw-r--r-- 1 user group 104745615 May  5 15:27 G1E_ATAC_0h_rep2.chr11.bam
-rw-r--r-- 1 user group 104157142 May  5 15:27 G1E_ATAC_30h_rep1.chr11.bam
-rw-r--r-- 1 user group   9254482 May  5 15:27 G1E_ATAC_30h_rep2.chr11.bam
-rw-r--r-- 1 user group  48952433 May  5 15:27 G1E_ChIP_Gata1_24h.chr11.bam
-rw-r--r-- 1 user group  20428179 May  5 15:27 G1E_ChIP_Gata2_0h.chr11.bam
-rw-r--r-- 1 user group  47261386 May  5 15:27 G1E_ChIP_Input_0h.chr11.bam
-rw-r--r-- 1 user group  37927946 May  5 15:27 G1E_ChIP_Input_24h.chr11.bam
~~~
{: .output}

We now need to make a new file that tells us how many reads there are mapping
to any given genomic position. One very common type of such file is a BigWig.
We can make these using a set of software called [deeptools][deeptools], which
we need to load on the cluster just as we did for bowtie2:

~~~
module avail
module load deeptools/3.0.1
bamCoverage --help
~~~
{: .language-bash}

We can get a detailed description of the options on the [bamCoverage help page][bamcoverage].

The command we'll want to use to map our ChIP and ATAC-seq files is:

~~~
bamCoverage -b aligned_reads/G1E_ATAC_0h_rep1.bam --binSize 1 --blackListFileName annotation/mm10-blacklist.ENCSR636HFF.v2.bed --normalizeUsing RPKM --extendReads -o bigwigs/G1E_ATAC_0h_rep1.bw
~~~
{: .language-bash}

~~~
'aligned_reads/G1E_ATAC_0h_rep1.bam' does not appear to have an index. You MUST index the file first!
~~~
{: .error}

Deeptools requires an index file. We can generate one of these with `samtools index`

~~~
samtools index aligned_reads/G1E_ATAC_0h_rep1.bam
~~~
{: .language-bash}

~~~
[E::hts_idx_push] unsorted positions
samtools index: "aligned_reads/G1E_ATAC_0h_rep1.bam" is corrupted or unsorted
~~~
{: .error}

We get another error, this time because `samtools index` only works on a file where
the reads are sorted by their chromosomal position. We can sort the files using
`samtools sort`

~~~
samtools sort aligned_reads/G1E_ATAC_0h_rep1.bam > aligned_reads/G1E_ATAC_0h_rep1.sorted.bam
samtools index aligned_reads/G1E_ATAC_0h_rep1.sorted.bam
ls -lh aligned_reads/
~~~
{: .language-bash}

~~~
-rw-r--r-- 1 user group  81M May  4 17:06 G1E_ATAC_0h_rep1.bam
-rw-r--r-- 1 user group  75M May  5 10:15 G1E_ATAC_0h_rep1.sorted.bam
-rw-r--r-- 1 user group  95K May  5 10:17 G1E_ATAC_0h_rep1.sorted.bam.bai
-rw-r--r-- 1 user group 100M May  4 16:58 G1E_ATAC_0h_rep2.bam
-rw-r--r-- 1 user group 100M May  4 17:11 G1E_ATAC_30h_rep1.bam
-rw-r--r-- 1 user group 8.9M May  4 17:12 G1E_ATAC_30h_rep2.bam
-rw-r--r-- 1 user group  47M May  4 16:43 G1E_ChIP_Gata1_24h.chr11.bam
-rw-r--r-- 1 user group  20M May  4 16:44 G1E_ChIP_Gata2_0h.chr11.bam
-rw-r--r-- 1 user group  46M May  4 16:45 G1E_ChIP_Input_0h.chr11.bam
-rw-r--r-- 1 user group  37M May  4 16:46 G1E_ChIP_Input_24h.chr11.bam
~~~
{: .output}

Now that we have sorted and indexed our bam file, we can finally run
bamCoverage:

~~~
bamCoverage -b aligned_reads/G1E_ATAC_0h_rep1.sorted.bam --binSize 1 \
--blackListFileName annotation/mm10-blacklist.ENCSR636HFF.v2.bed \
--normalizeUsing RPKM --extendReads -o bigwigs/G1E_ATAC_0h_rep1.bw -r chr11:32180000-32345000
ls -l bigwigs/
~~~
{: .language-bash}

> ## Sorting and indexing files in a loop (part 1)
>
> Write a pattern to match just the unsorted .bam files
>
> ~~~
> ls -l _______
> ~~~
> {: .language-bash}
>
> ### Expected output:
>
> ~~~
> -rw-r--r-- 1 user group  84893348 May  4 17:06 aligned_reads/G1E_ATAC_0h_rep1.chr11.bam
> -rw-r--r-- 1 user group 104745615 May  4 16:58 aligned_reads/G1E_ATAC_0h_rep2.chr11.bam
> -rw-r--r-- 1 user group 104157142 May  4 17:11 aligned_reads/G1E_ATAC_30h_rep1.chr11.bam
> -rw-r--r-- 1 user group   9254482 May  4 17:12 aligned_reads/G1E_ATAC_30h_rep2.chr11.bam
> -rw-r--r-- 1 user group  48952433 May  4 16:43 aligned_reads/G1E_ChIP_Gata1_24h.chr11.bam
> -rw-r--r-- 1 user group  20428179 May  4 16:44 aligned_reads/G1E_ChIP_Gata2_0h.chr11.bam
> -rw-r--r-- 1 user group  47261386 May  4 16:45 aligned_reads/G1E_ChIP_Input_0h.chr11.bam
> -rw-r--r-- 1 user group  37927946 May  4 16:46 aligned_reads/G1E_ChIP_Input_24h.chr11.bam
> ~~~
> {: .output}
>
> > ## Solution
> >
> > ~~~
> > ls -l aligned_reads/G1E_*.chr11.bam
> > ~~~
> > {: .language-bash}
> {: .solution}
{: .challenge}

> ## Sorting and indexing files in a loop (part 2)
>
> Use the pattern you just wrote to write a loop that uses the
> `cut` command to print out the name of the bam file without the
> ".bam" part.
>
> ~~~
> for bamfile in <PATTERN FROM PART 1>
> do
>   echo $bamfile
>   echo $bamfile | cut -d ____ -f ____
> done
> ~~~
> {: .language-bash}
>
> ### Expected output:
>
> ~~~
> aligned_reads/G1E_ATAC_0h_rep1.chr11.bam
> aligned_reads/G1E_ATAC_0h_rep1.chr11
> aligned_reads/G1E_ATAC_0h_rep2.chr11.bam
> aligned_reads/G1E_ATAC_0h_rep2.chr11
> aligned_reads/G1E_ATAC_30h_rep1.chr11.bam
> aligned_reads/G1E_ATAC_30h_rep1.chr11
> aligned_reads/G1E_ATAC_30h_rep2.chr11.bam
> aligned_reads/G1E_ATAC_30h_rep2.chr11
> aligned_reads/G1E_ChIP_Gata1_24h.chr11.bam
> aligned_reads/G1E_ChIP_Gata1_24h.chr11
> aligned_reads/G1E_ChIP_Gata2_0h.chr11.bam
> aligned_reads/G1E_ChIP_Gata2_0h.chr11
> aligned_reads/G1E_ChIP_Input_0h.chr11.bam
> aligned_reads/G1E_ChIP_Input_0h.chr11
> aligned_reads/G1E_ChIP_Input_24h.chr11.bam
> aligned_reads/G1E_ChIP_Input_24h.chr11
> ~~~
> {: .output}
> > ## Solution
> >
> > ~~~
> for bamfile in aligned_reads/G1E_*.chr11.bam
> do
>   echo $bamfile
>   echo $bamfile | cut -d "." -f 1,2
> done
> > ~~~
> > {: .language-bash}
> {: .solution}
{: .challenge}

> ## Sorting and indexing files in a loop (part 3)
>
> Using your answers to parts 1 & 2 and the "command substitution"
> syntax `$()` we learned in the [shell find lesson][shell-find]
> to print the original bam file and its corresponding
> sorted bam file.
>
> ~~~
> for bamfile in <PATTERN FROM PART 1>
> do
>   echo $bamfile $(echo $bamfile | cut -d ____ -f ____).sorted.bam
> done
> ~~~
> {: .language-bash}
>
> ### Expected output:
>
> ~~~
> aligned_reads/G1E_ATAC_0h_rep1.chr11.bam aligned_reads/G1E_ATAC_0h_rep1.chr11.sorted.bam
> aligned_reads/G1E_ATAC_0h_rep2.chr11.bam aligned_reads/G1E_ATAC_0h_rep2.chr11.sorted.bam
> aligned_reads/G1E_ATAC_30h_rep1.chr11.bam aligned_reads/G1E_ATAC_30h_rep1.chr11.sorted.bam
> aligned_reads/G1E_ATAC_30h_rep2.chr11.bam aligned_reads/G1E_ATAC_30h_rep2.chr11.sorted.bam
> aligned_reads/G1E_ChIP_Gata1_24h.chr11.bam aligned_reads/G1E_ChIP_Gata1_24h.chr11.sorted.bam
> aligned_reads/G1E_ChIP_Gata2_0h.chr11.bam aligned_reads/G1E_ChIP_Gata2_0h.chr11.sorted.bam
> aligned_reads/G1E_ChIP_Input_0h.chr11.bam aligned_reads/G1E_ChIP_Input_0h.chr11.sorted.bam
> aligned_reads/G1E_ChIP_Input_24h.chr11.bam aligned_reads/G1E_ChIP_Input_24h.chr11.sorted.bam
> ~~~
> {: .output}
> > ## Solution
> >
> > ~~~
> for bamfile in aligned_reads/G1E_*.chr11.bam
> do
>   echo $bamfile $(echo $bamfile | cut -d "." -f 1,2).sorted.bam
> done
> > ~~~
> > {: .language-bash}
> {: .solution}
{: .challenge}

> ## Sorting and indexing files in a loop (part 4)
>
> Using your answers to parts 1, 2 & 3 write a for loop that sorts
> and indexes all of the bam files.
> ~~~
> for ______________
> do
>   samtools sort ___________
>   samtools index __________
> done
> ~~~
> {: .language-bash}
>
> > ## Solution
> >
> > ~~~
> for bamfile in aligned_reads/G1E_*.chr11.bam
> do
>   samtools sort $bamfile > $(echo $bamfile | cut -d "." -f 1,2).sorted.bam
>   samtools index $(echo $bamfile | cut -d "." -f 1,2).sorted.bam
> done
> > ~~~
> > {: .language-bash}
> {: .solution}
{: .challenge}

And since we've now sorted and indexed all our bam files, we should be able
to make bigwigs for them all using a for loop:

~~~
for bam in aligned_reads/G1E_*.sorted.bam
do
  bamCoverage -b $bam --binSize 1 --blackListFileName annotation/mm10-blacklist.ENCSR636HFF.v2.bed \
  --normalizeUsing RPKM --extendReads -o $(echo $bam | cut -d "." -f 1,2).bw -r chr11:32180000-32345000
  mv $(echo $bam | cut -d "." -f 1,2).bw bigwigs/
done
~~~
{: .language-bash}

~~~
*ERROR*: library is not paired-end. Please provide an extension length.
mv: cannot stat `aligned_reads/G1E_ChIP_Gata2_0h.chr11.bw': No such file or directory
~~~
{: .error}

We get an error message because the ChIP samples are not paired end, so we can't
automatically extend the fragments, we have to give an estimated fragment size:

~~~
for bam in aligned_reads/G1E_ChIP*.sorted.bam
do
  bamCoverage -b $bam --binSize 1 --blackListFileName annotation/mm10-blacklist.ENCSR636HFF.v2.bed \
  --normalizeUsing RPKM --extendReads 150 -o $(echo $bam | cut -d "." -f 1,2).bw -r chr11:32180000-32345000
  mv $(echo $bam | cut -d "." -f 1,2).bw bigwigs/
done
ls -l bigwigs/
~~~
{: .language-bash}

~~~
-rw-r--r-- 1 user group 24116 May  5 15:01 G1E_ATAC_0h_rep1.chr11.bw
-rw-r--r-- 1 user group 28656 May  5 15:02 G1E_ATAC_0h_rep2.chr11.bw
-rw-r--r-- 1 user group 28260 May  5 15:02 G1E_ATAC_30h_rep1.chr11.bw
-rw-r--r-- 1 user group  4031 May  5 15:02 G1E_ATAC_30h_rep2.chr11.bw
-rw-r--r-- 1 user group 38176 May  5 15:07 G1E_ChIP_Gata1_24h.chr11.bw
-rw-r--r-- 1 user group 24922 May  5 15:07 G1E_ChIP_Gata2_0h.chr11.bw
-rw-r--r-- 1 user group 33703 May  5 15:07 G1E_ChIP_Input_0h.chr11.bw
-rw-r--r-- 1 user group 27264 May  5 15:08 G1E_ChIP_Input_24h.chr11.bw
~~~
{: .output}

{% include links.md %}

