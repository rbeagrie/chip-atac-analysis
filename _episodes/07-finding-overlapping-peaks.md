---
title: "Finding overlapping peaks"
teaching: 0
exercises: 0
questions:
- "Key question (FIXME)"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

For this session, we're going to be using `bedtools`. We need to load it in to our shell first:

~~~
module load bedtools/2.27.1
bedtools intersect --help
~~~
{: .language-bash}

You should use the peak files ending in ".narrowPeak" for the following challenges.

> ## Getting peak calls
>
> If you didn't manage to finish peak calling all of the files from lesson 6, you can download some
> pre-computed peak calls and use those instead
>
> ~~~
> cd chip-tutorial
> wget https://rob.beagrie.com/media/chip-tutorial/peak_calls.tar.gz
> tar zxvf peak_calls.tar.gz
> ls peak_calls/macs2/*.narrowPeak
> ~~~
> {: .source}
>
> ~~~
> peak_calls/macs2/G1E_ATAC_0h_rep1_peaks.narrowPeak
> peak_calls/macs2/G1E_ATAC_0h_rep2_peaks.narrowPeak
> peak_calls/macs2/G1E_ATAC_30h_rep1_peaks.narrowPeak
> peak_calls/macs2/G1E_ATAC_30h_rep2_peaks.narrowPeak
> peak_calls/macs2/G1E_ChIP_Gata1_24h_peaks.narrowPeak
> peak_calls/macs2/G1E_ChIP_Gata2_0h_peaks.narrowPeak
> ~~~
> {: .output}
{: .callout}

> ## Challenge 1: How many overlapping peaks are there?
>
> Use `bedtools intersect` to work out how many peaks overlap between replicates 1 and 2 of the
> ATAC-seq data from 0h (uninduced) G1E cells. You can use `bedtools intersect --help`, the
> [bedtools online help page][bedtools-intersect] and google for help.
>
> > ## Solution
> >
> > ~~~
> > bedtools intersect -u -a peak_calls/macs2/G1E_ATAC_0h_rep1_peaks.narrowPeak \
> > -b peak_calls/macs2/G1E_ATAC_0h_rep2_peaks.narrowPeak | wc -l
> > ~~~
> > {: .language-bash}
> >
> {: .solution}
{: .challenge}


> ## Challenge 2: How many overlapping peaks are there (part 2)?
>
> Write a nested for loop (i.e. two for loops, one within the other) to count the number of
> overlapping peaks between every possible combination of the peak files from the ATAC-seq
> samples.
>
> HINT: Remember to build your loop commands up slowly, starting from the smallest possible piece of
> code and slowly adding things on.
>
> > ## Solution
> >
> > ~~~
> > for bed1 in peak_calls/macs2/G1E_ATAC_*.narrowPeak
> > do
> >   for bed2 in peak_calls/macs2/G1E_ATAC_*.narrowPeak
> >   do
> >     echo $(basename $bed1) $(basename $bed2) $(bedtools intersect -u -a $bed1 -b $bed2 | wc -l)
> >   done
> > done
> > ~~~
> > {: .language-bash}
> >
> {: .solution}
{: .challenge}


> ## Challenge 3: How many overlapping peaks are there (part 3)?
>
> Hopefully you've found out that there is quite a lot of overlap between our peak calls. But how
> many distinct peaks did we call in total? We can find out using [bedtools merge][bedtools-merge].
> You will also need to use `bedtools sort`
>
> HINT: It might be easier to make an intermediate file at each step, but you can also do this as a
> single command in a pipeline.
>
> > ## Solution
> >
> > ~~~
> > cat peak_calls/macs2/G1E_ATAC_*.narrowPeak | bedtools sort | bedtools merge | wc -l
> > ~~~
> > {: .language-bash}
> >
> {: .solution}
{: .challenge}

> ## Challenge 4: Making a list of consensus ATAC-seq peaks
>
> We know how many distinct peaks we had in total, but we probably don't want to trust those that
> were only identified in a single replicate. We need to get bedtools merge to retain information
> about which peaks were merged together. If you look in one of the narrowPeak files you will
> see that column 4 is potentially helpful here:
>
> ~~~
> head peak_calls/macs2/G1E_ATAC_0h_rep1_peaks.narrowPeak
> ~~~
> {: .language-bash}
>
> ~~~
> chr11	3104874	3105821	G1E_ATAC_0h_rep1_peak_1	291	.	9.47788	32.42815	29.18559	557
> chr11	3112182	3112345	G1E_ATAC_0h_rep1_peak_2	12	.	3.02204	3.25904	1.28949	109
> chr11	3123275	3124132	G1E_ATAC_0h_rep1_peak_3	829	.	15.12556	89.32738	82.92233	600
> chr11	3138607	3138829	G1E_ATAC_0h_rep1_peak_4	210	.	13.28701	24.07859	21.09448	101
> chr11	3152358	3152962	G1E_ATAC_0h_rep1_peak_5	119	.	3.50204	14.67430	11.98344	234
> chr11	3153981	3154339	G1E_ATAC_0h_rep1_peak_6	51	.	2.45218	7.50518	5.10658	125
> chr11	3155042	3156271	G1E_ATAC_0h_rep1_peak_7	63	.	2.62152	8.83316	6.37181	1001
> chr11	3156443	3156962	G1E_ATAC_0h_rep1_peak_8	127	.	3.28824	15.49932	12.78019	131
> chr11	3157227	3157683	G1E_ATAC_0h_rep1_peak_9	55	.	2.46292	7.97198	5.55084	26
> chr11	3157856	3158101	G1E_ATAC_0h_rep1_peak_10	41	.	2.35471	6.49318	4.14949	90
> ~~~
> {: .output}
>
> Read the [bedtools merge][bedtools-merge] documentation and figure out how to report all the
> distinct values of column 4 from the peaks that were merged, and the number of distinct values
> that were merged ('distinct' values and 'count_distinct'). Make sure to save the results in a file
> with the merged peak calls.
>
> > ## Solution
> >
> > ~~~
> > cat peak_calls/macs2/G1E_ATAC_*.narrowPeak | bedtools sort | bedtools merge -c 4,4 -o \
> > distinct,count_distinct > peak_calls/macs2/G1E_ATAC_all_merged.bed
> > ~~~
> > {: .language-bash}
> >
> {: .solution}
>
> **BONUS:** the command `cat peak_calls/macs2/G1E_ATAC_*.narrowPeak | sed 's/_peak_/\t/'` will split
> column 4 into two colums, so that the new column 4 contains only the name of the sample the peak was
> found in and the new column 5 contains a unique identifier for each peak from that file. Why would
> this be helpful to add to your pipeline?
>
{: .challenge}

> ## Challenge 5: Making a list of consensus ATAC-seq peaks (part 2)
>
> As long as the number of peaks that were merged together is the last column in your file,
> `grep "1$"` will give you a list of peaks that were only found in one replicate (this works
> because the $ symbol indicates the end of a line, so grep will only find lines that end with a 1
> as the last character). How can you modify this to find peaks that were found in more than one
> replicate and save them to a new file?
>
> > ## Solution
> >
> > ~~~
> > grep -v "1$" peak_calls/macs2/G1E_ATAC_all_merged.bed > \
> > peak_calls/macs2/G1E_ATAC_all_reproducible.bed
> > ~~~
> > {: .language-bash}
> >
> {: .solution}
>
> **BONUS:** This command wouldn't work if you were merging peaks from more than ten input files. Why?
> It also doesn't allow you to increase the threshold of the minimum number of samples a peak needs
> to be found in to be considered a real peak. One alternative approach would be to import the peaks
> file as a pandas Dataframe using Python, and to write a python script that would filter peaks in a
> more sophisticated manner. You can have a go at this if you have some extra time.
>
{: .challenge}

> ## Challenge 6: How big are our peaks?
>
> Now that we have a list of consensus peaks, we can start thinking about doing an analysis to look
> for differential accessibility. In order to do this, we need a table telling us how many reads
> mapped to each peak (rows) in each sample (columns). Read the help page for [bedtools
> multicov][bedtools-multicov] to work out how to generate this table.
>
> > ## Solution
> >
> > ~~~
> > bedtools multicov -bed peak_calls/macs2/G1E_ATAC_all_reproducible.bed \
> > -bams aligned_reads/G1E_ATAC_0h_rep1.chr11.sorted.bam \
> > aligned_reads/G1E_ATAC_0h_rep2.chr11.sorted.bam \
> > aligned_reads/G1E_ATAC_30h_rep1.chr11.sorted.bam \
> > aligned_reads/G1E_ATAC_30h_rep2.chr11.sorted.bam > \
> > peak_counts_by_sample.table
> > ~~~
> > {: .language-bash}
> >
> > We don't want to write `-bams aligned_reads/G1E_ATAC_*.sorted.bam` because then we won't know
> > what order the columns were in for downstream analysis. In fact it would be even better if this
> > new table had a header line so we can remember which column was which:
> >
> > ~~~
> > bams=$(echo aligned_reads/G1E_ATAC_*.sorted.bam)
> > echo "chrom  start  stop  $bams" > peak_counts_by_sample.table
> > bedtools multicov -bed peak_calls/macs2/G1E_ATAC_all_reproducible.bed \
> > -bams $bams >> peak_counts_by_sample.table
> > ~~~
> > {: .language-bash}
> >
> >
> {: .solution}
>
> **BONUS:** Now that we have a table of peak heights for each sample, we can see how similar our
> biological replicates are. Open the table in Python as a pandas dataframe and make scatter plots
> to compare the samples. Are biological replicates more similar to one another than to the other
> time point (i.e. is 0h rep1 more similar to 0h rep2 or 30h rep1)?
>
{: .challenge}

{% include links.md %}

