---
title: "Motif analysis"
teaching: 0
exercises: 0
questions:
- "Key question (FIXME)"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

### Looking for occurrences of known motifs

In this session, we will examine how to find motifs in our called peaks.

The first thing we need to do is to upload the peaks that we called with DESeq2 back to the cluster.

Do this using scp from your local machine:

~~~
scp Downloads/peaks_significantly_down_at_30h.bed user@cluster:~/chip-tutorial/
scp Downloads/peaks_significantly_up_at_30h.bed user@cluster:~/chip-tutorial/
~~~
{: .language-bash}

In order to find motifs in these peaks, we first need their sequences. We can fetch these using
`bedtools `:

~~~
cd chip-tutorial
mkdir motif_analysis
module load bedtools/2.27.1
bedtools getfasta -fi /path/to/fasta_file -bed peaks_significantly_up_at_30h.bed > \
motif_analysis/peaks_significantly_up_at_30h.fa
head -n 2 motif_analysis/peaks_significantly_up_at_30h.fa
~~~
{: .language-bash}

~~~
>chr11:3648419-3649280
ggcccggaagcaaacccagatcttccgactccGGAGCTCGAATCCCCACCCCGGAGCTGCCGCCAGCGCTCTCGCTCCCCACTAAGGTGGCATAAGTAGGAAAAGGCAAAGGAGGAAGCCCACCACCAACCATCCTCCAGCTCCCTCTTGACCCTCCAACTGGGTCCTAGCAGGCCCCGAGGGGTGGCGCCAGGCCTCCAGGCCAGCCCCGGCCACCCCAGCTTCCCCCGCCGCCTCCCCGGGAGCGGGGCGGAGTTGCCCGCCGCCGCCAACCACCCGGACGACCCCTCGTCGGATCGCAAAAGCATAGGCCGCCGCCCGAGTTCTGCGTACGAGAAGAAAGACGCGGCGCGAGCGCCAACGGCCAccgggcgcgcgccgcggcggccgggcctgcgccccAAGAGCTGGATGCCAGAGCAGGAGAAAGAGCCGCCAACCGATCGCTCGATCGACCGCACGCCCGCTCCTTCGCCCTACCAGACCGGGGAGGGGGGGAGGGCGCGCCAGGGCTTCTTCGAGTTAGGGGCCTGACTCCCCGGCGACAAGACAAGATGGCTCCTCCGCTCTGCTCTGCAGCTACCGTCCCGGCTGCGTCTCTCGGCCCGCCCCCGGTACTGTTCCTTTAAATGGCCGAGAGGGGTCAGGAAAGCAGGAAGAGAGGCGCGCGCCCTCACGTGATGGGGGCGGACCGGCCAAACAAGATCCGGTCACGTGGACCGACGTCGCTCGCCAATCAGAAGGTGGAGTCTCAGGGCGGGGGCGGGAAGACTTCGAGCTCTGTGATTGGCCGAAGGTGCGCTGGAGGCCTCGGCGGCCGGGCCGTCGGGAGTGAGCGAGGTGATTCGCGCTGCGTGCGCT
~~~
{: .output}

> ## Challenge: Fetching peak sequences
>
> Now make a file containing the sequences of the peaks that had lower accessibility at 30h.
>
> > ## Solution
> >
> > ~~~
> > bedtools getfasta -fi /path/to/fasta_file -bed peaks_significantly_down_at_30h.bed > \
> > motif_analysis/peaks_significantly_down_at_30h.fa
> > ~~~
> > {: .language-bash}
> >
> {: .solution}
>
{: .challenge}

We have the peak sequences, so now we need a motif to search for. Go to [the hocomoco motif
database][hocomoco] and search for **Gata1** (make sure to select the mouse database). We're
presented with two options, a short core motif and a longer motif with the Tal1 upstream site, so
click on the shorter motif, then scroll down to find the download link for the *pcm* file and right
click to copy the link address.

Once you have the address, you can download the motif file directly to the cluster using `wget`:

~~~
cd motif_analysis
wget https://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/MOUSE/mono/pcm/GATA1_MOUSE.H11MO.1.A.pcm
~~~
{: .language-bash}

We're going to use `fimo` to search for Gata1 motif occurrences in our peaks, but fimo requires a
slightly different format for the motif file. Luckily it comes with a conversion utility we can use
to create the right sort of motif file. Fimo is part of the meme suite, so we need to find an
appropriate module to load (for some reason meme version 5 and above don't seem to come with motif
conversion tools, so look for a module for meme v4).

~~~
module avail
module load meme/4.9.1_1
chen2meme GATA1_MOUSE.H11MO.1.A.pcm > GATA1_MOUSE.H11MO.1.A.meme
~~~
{: .language-bash}

Now that we have a motif file, we can use it with fimo:

~~~
fimo -o gata1_up_peaks GATA1_MOUSE.H11MO.1.A.meme peaks_significantly_up_at_30h.fa
ls -l gata1_up_peaks
~~~
{: .language-bash}

~~~
-rw-r--r-- 1 user group   265 May 19 09:13 cisml.css
-rw-r--r-- 1 user group 32511 May 19 09:13 cisml.xml
-rw-r--r-- 1 user group  6531 May 19 09:13 fimo.gff
-rw-r--r-- 1 user group 19520 May 19 09:13 fimo.html
-rw-r--r-- 1 user group 14261 May 19 09:13 fimo-to-html.xsl
-rw-r--r-- 1 user group  3634 May 19 09:13 fimo.txt
-rw-r--r-- 1 user group  1598 May 19 09:13 fimo.xml
~~~
{: .output}

You can have a look at these files to see if you can figure out what information they contain. For
now, let's examine gata1_up_peaks/fimo.gff:

~~~
head -n5 gata1_up_peaks/fimo.gff
~~~
{: .language-bash}

~~~
##gff-version 3
chr11:101732559-101733190	fimo	nucleotide_motif	58	68	51.0	+	.	Name=GATA1_MOUSE.H11MO.1.A;ID=GATA1_MOUSE.H11MO.1.A-1;pvalue=8e-06;qvalue=0.366;sequence=ACAGATAAGCA
chr11:102363673-102363981	fimo	nucleotide_motif	82	92	42.6	+	.	Name=GATA1_MOUSE.H11MO.1.A;ID=GATA1_MOUSE.H11MO.1.A-2;pvalue=5.53e-05;qvalue=0.704;sequence=CCTGATAAGAG
chr11:102363673-102363981	fimo	nucleotide_motif	167	177	56.6	+	.	Name=GATA1_MOUSE.H11MO.1.A;ID=GATA1_MOUSE.H11MO.1.A-1;pvalue=2.17e-06;qvalue=0.198;sequence=ACAGATAAGGG
chr11:103412015-103413050	fimo	nucleotide_motif	250	260	49.3	+	.	Name=GATA1_MOUSE.H11MO.1.A;ID=GATA1_MOUSE.H11MO.1.A-1;pvalue=1.18e-05;qvalue=0.394;sequence=TGAGATAAGAG
~~~
{: .output}

The first line is a header, then each subsequent line of this gff file is a match to a Gata1 motif that fimo found.

> ## Challenge: Counting motif occurrances
>
> Use fimo to look for Gata1 motifs in the peaks with decreased accessibility. Is there a difference
> in the frequency of Gata1 motifs between the "up" peaks and the "down" peaks?
> or the "down" peaks
>
> > ## Solution
> >
> > ~~~
> > fimo -o gata1_down_peaks GATA1_MOUSE.H11MO.1.A.meme peaks_significantly_down.fa
> > wc -l gata1_*_peaks/fimo.gff
> > ~~~
> > {: .language-bash}
> >
> > ~~~
> >   190 gata1_down_peaks/fimo.gff
> >   40 gata1_up_peaks/fimo.gff
> >  230 total
> > ~~~
> > {: .output}
> >
> > Interestingly, there are a lot more Gata1 motifs in the peaks that decrease in accessibility
> > than those which increase in accessibility, even though there are similar numbers of peaks in
> > each list.
> >
> {: .solution}
>
{: .challenge}

### De novo motif analysis

If we don't know what motifs we want to look for, we can do a *de novo* motif analysis that tries to
generate motifs which are enriched in a set of peaks and then match those back to known
transcription factors. The package we use for this is called [meme chip][meme-chip], which we've
already loaded as part of the meme suite.

~~~
qlogin
module load meme/5.0.1
meme-chip -o motif_analysis/meme_up_peaks -db /path/to/meme_database \
motif_analysis/peaks_significantly_up_at_30h.fa
~~~
{: .language-bash}

To look through the results, we'll need to download the whole folder to our local machine. We can do
that with `scp -r`:

~~~
scp -r user@cluster:/home/user/chip-tutorial/motif_analysis/meme_up_peaks .
cd meme_up_peaks
open meme.html
~~~
{: .language-bash}



{% include links.md %}

