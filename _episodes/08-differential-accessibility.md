---
title: "Calling differential accessibility using DESeq2"
teaching: 0
exercises: 0
questions:
- "Key question (FIXME)"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

We're going to use DESeq2 on our local machines, so first we need to download the peaks file from
the previous lesson using `scp` and the absolute path of the file.

~~~
ls $(pwd)/peak_counts_by_sample.table
~~~
{: .language-bash}

~~~
/home/user/chip-tutorial/peak_counts_by_sample.table
~~~
{: .output}

~~~
scp user@cluster:/home/user/chip-tutorial/peak_counts_by_sample.table .
~~~
{: .language-bash}


~~~
mkdir peak_calls
mkdir peak_calls/macs2
macs2 callpeak -t aligned_reads/G1E_ChIP_Gata2_0h.chr11.sorted.bam \
-c aligned_reads/G1E_ChIP_Input_0h.chr11.sorted.bam -n G1E_ChIP_Gata2_0h \
-f BAM -g 1.87e9 -q 0.1 --outdir peak_calls/macs2
ls -l peak_calls/macs2
~~~
{: .language-bash}

Once macs2 has finished running, we should have four new files in our output directory.

~~~
-rw-r--r-- 1 user group  88633 May 12 12:06 G1E_ChIP_Gata2_0h_model.r
-rw-r--r-- 1 user group 135241 May 12 12:06 G1E_ChIP_Gata2_0h_peaks.narrowPeak
-rw-r--r-- 1 user group 151285 May 12 12:06 G1E_ChIP_Gata2_0h_peaks.xls
-rw-r--r-- 1 user group  93538 May 12 12:06 G1E_ChIP_Gata2_0h_summits.bed
~~~
{: .output}

We can see that the ".narrowPeak" file is in BED format with a few extra columns.
You can read about this file format at [UCSC][narrowpeak-format].
~~~
head peak_calls/macs2/G1E_ChIP_Gata2_0h_peaks.narrowPeak
~~~
{: .language-bash}

~~~
chr11	3104926	3105658	G1E_ChIP_Gata2_0h_peak_1	998	.	11.03180	103.14955	99.81662	416
chr11	3105932	3106226	G1E_ChIP_Gata2_0h_peak_2	61	.	3.02292	8.60469	6.17999	159
chr11	3201865	3202318	G1E_ChIP_Gata2_0h_peak_3	105	.	3.96565	13.13776	10.59181	374
chr11	3279535	3280209	G1E_ChIP_Gata2_0h_peak_4	628	.	16.61709	65.97417	62.87525	400
chr11	3291562	3291949	G1E_ChIP_Gata2_0h_peak_5	102	.	4.20686	12.83781	10.29883	159
chr11	3313151	3313768	G1E_ChIP_Gata2_0h_peak_6	1158	.	24.21939	119.31017	115.88039	324
chr11	3379576	3379870	G1E_ChIP_Gata2_0h_peak_7	55	.	5.12864	7.99220	5.58688	90
chr11	3434737	3435288	G1E_ChIP_Gata2_0h_peak_8	175	.	9.95468	20.25646	17.57889	272
chr11	3446439	3446733	G1E_ChIP_Gata2_0h_peak_9	46	.	4.69043	6.99010	4.62133	264
chr11	4086358	4086652	G1E_ChIP_Gata2_0h_peak_10	49	.	4.89368	7.31591	4.93340	69
~~~
{: .output}

~~~
macs2 callpeak -t aligned_reads/G1E_ChIP_Gata1_24h.chr11.sorted.bam \
-c aligned_reads/G1E_ChIP_Input_24h.chr11.sorted.bam -n G1E_ChIP_Gata1_24h \
-f BAM -g 1.87e9 -q 0.1 --outdir peak_calls/macs2
wc -l peak_calls/macs2/*.narrowPeak
~~~
{: .language-bash}

~~~
   9 peak_calls/macs2/G1E_ChIP_Gata1_24h_peaks.narrowPeak
1548 peak_calls/macs2/G1E_ChIP_Gata2_0h_peaks.narrowPeak
1557 total
~~~
{: .output}

The Gata1 file has only 9 peaks. We should examine them somehow to see what's going on.
We can upload them to the UCSC genome browser, but first we need to download them to our
local computer, which we can do using `scp`. We need the absolute path to the file
we need to download, which we can get using a combination of `ls` and `pwd`.

Once we have the absolute path, we can open another terminal and download the files
like this:

~~~
scp username@cluster.com:/home/user/chip-tutorial/peak_calls/macs2/*.narrowPeak Downloads/
~~~
{: .language-bash}

Then we copy and paste the data into [UCSC][ucsc-mm10]. Open the genome browser and go to "add
custom tracks". We need to paste our data into the text box, with a track line first.
For Gata1, this might look like:

~~~
track type=narrowPeak name="Gata1-ChIP" description="ChIP for Gata1 in 24h induced G1E cells"
~~~
{: .language-bash}

We can see the reason macs2 hasn't called many peaks for Gata1 is that it looks very similar
to the control, input track (i.e. it's not a very good ChIP).

Now we can peak call the ATAC-seq datasets. We need to make two changes to the command for
these samples because there's no control track and because they are paired end files.

~~~
macs2 callpeak -t aligned_reads/G1E_ATAC_0h_rep1.chr11.sorted.bam \
-n G1E_ATAC_0h_rep1 -f BAMPE -g 1.87e9 -q 0.1 --outdir peak_calls/macs2
~~~
{: .language-bash}

Let's write a batch script to peak call all of the ATAC-seq files. First log out of the compute
node using <kbd>Ctrl</kbd>+<kbd>D</kbd>

> ## Peak calling the rest of the ATAC-seq files.
>
> Open a new bash script called `peak_call_atac.sh` and use the template below to write a script for
> submission to the cluster that will peak call all of the ATAC-seq files.
>
> ~~~
> #$ -cwd
> #$ -m ea
>
> for ____
> do
>   macs2 callpeak ____ -n $(basename $bam | cut -d "." -f 1)
> done
> ~~~
> {: .source}
>
> > ## Solution
> >
> > ~~~
> > #$ -cwd
> > #$ -m ea
> >
> > for bam in aligned_reads/G1E_ATAC_*.sorted.bam
> > do
> >  macs2 callpeak -t $bam -n $(basename $bam | cut -d "." -f 1) \
> >  -f BAMPE -g 1.87e9 -q 0.1 --outdir peak_calls/macs2
> > done
> > ~~~
> > {: .source}
> >
> {: .solution}
{: .challenge}

{% include links.md %}

