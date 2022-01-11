# TimemaFusion

Characterizing chromosomal fusion(s) in *Timema* stick insects

## Data sets

Most of the data sets for this project are in `/uufs/chpc.utah.edu/common/home/gompert-group3/data/timema_clines_rw_SV/`; the nanopore data for *T. cristinae* from Refugio are in `/uufs/chpc.utah.edu/common/home/gompert-group3/data/nanoporereads/TimemaL1`.

0. Whole genomes for compariative genome alignment, HiC genomes = *T. cristinae* striped, *T. cristinae* green unstriped, *T. knulli*, *T. podura* and *T. chumash* (don't necessarily care about last two for this project).
1. **hwy154** = GBS data, LD for comparison with Refugio to test of fusion
2. **fha_mapping_sample** = GBS data, LD as well?
3. **n1_doro** = GBS data, LD as well?
4. **refugio** and **refugio_2019** = GBS data, LD to test for fusion
5. **rw_plus** = GBS data, plans for mapping performance, patterns of LD, ect.
6. **refugio_nanoper** = nanpore data from 6 *Timema cristinae* from Refugio, could provide direct evidence of a fusion on Refugio. 
7. **matepair_SVs** = matepair data form Kay's paper, could provide direct evidence of a fusion on Refugio.

## Comparative genome alignments

From past alignments (details to come) Green striped 12033 is LG11 (which we are thinking about for the redwood stuff) and Green 7748 = LG8.

As a first, more general pass, I am assessing synteny (not colinearity) simply using blast. This seems to be working well.

* First, LGs from melanic *T. cristinae* versus the green striped *T. cristinae* HiC genome.

The 13 (unambiguous) large scaffold (> 10 Mbps) from the green striped genome are listed, along with their sizes in bps, in `/uufs/chpc.utah.edu/common/home/gompert-group1/data/timema/hic_genomes/t_crist/greenGenomeLGscafs.txt`. I used this to extract these scaffolds from the genome with `samtools` (version 1.12)

```{bash}
samtools faidx timema_cristinae_12Jun2019_lu3Hs.fasta -r greenLGRegions.txt -o timema_cristinae_LGs_12Jun2019_lu3Hs.fasta
```
Next I created a blastable data base to match scaffolds from the melanic genome to the green striped genome. This was done with `blast` (version 2.11.0).

```{bash}
makeblastdb -in timema_cristinae_LGs_12Jun2019_lu3Hs.fasta -dbtype nucl -parse_seqids -out GreenGenome -title "Tcr Green genome chrom. scafs."
```

Trying to use `blastn` to match chromosomes. First, I ran a strict blast search.

```{bash}
blastn -db GreenGenome -evalue 1e-50 -perc_identity 92 -query ../../tcrDovetail/version3/mod_map_timema_06Jun2016_RvNkF702.fasta -outfmt 6 -num_threads 48 > Green2Melanic.txt
```
Then, I worte a `perl` script to sum the total alignment length between each of the 13 scaffolds from the green striped genome and the 13 LGs we have defined for the melanic genome. I could also do this at the scaffold level (as it is unlikely that all scaffold assignments to LGs from the mapping families are correct), but this seems more useful as the main goal is to preserve as many of the LG IDs for the newer genome as possible to maximize consistency across papers (i.e., using the best overall match between old LGs and new scaffolds is sufficient). It is important (empirically) to only count long alignments; I tried counting all and you just get noise.

```{perl}
#!/usr/bin/perl
#
# compute total alignment lengths between melanic LGs and green stripe scaffolds
#

$blast = "Green2Melanic.txt";

## column headers, see blastn -help for details
## qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore


open(IN, $blast) or die "failed to read $blast\n";
open(OUT, "> AlnGreenStr2Melanic.txt") or die "failed to write OUT\n";

while(<IN>){
        chomp;
        $_ =~ m/LG\-(\d+)\S+\s+Sclu3Hs_(\d+)\S+\s+\S+\s+(\d+)/ or print "No match here: $_\n";
        $lg = $1;
        $sc = $2;
        $len = $3;
        $aln = "$sc"."_$lg";
        if($len > 10000){ ## only consider long alignments
                if(defined $alns{$aln}){
                        $alns{$aln} += $len;
                }
                else {
                        $alns{$aln} = $len;
                }
        }
}
close(IN);
foreach $aln (sort keys %alns){
        print OUT "$aln $alns{$aln}\n";
}
```

Last, I used `R` to summarize the alignments and map the melanic LGs to the big scaffolds (chromosomes) from the green striped *T. cristinae* genome.

```{R}
aln<-read.table("AlnGreenStr2Melanic.txt",header=FALSE)

amat<-matrix(as.numeric(unlist(strsplit(x=as.character(aln[,1]),"_"))),nrow=59,ncol=2,byrow=TRUE)

aln_mat<-matrix(0,nrow=13,ncol=14)
usc<-unique(amat[,1])
for(i in 1:13){for(j in 1:14){
    k<-j-1
    a<-which(amat[,1]==usc[i] & amat[,2]==k)
    if(length(a) == 1){
            aln_mat[i,j]<-aln[a,2]
    }
}}

colnames(aln_mat)<-0:13
row.names(aln_mat)<-usc

for(i in 1:13){
    aln_mat[i,]<-aln_mat[i,]/sum(aln_mat[i,])
    }


library(fields)
pdf("GreenStr2Melanic.pdf",width=6,height=6)
cs<-hcl.colors(10, "YlOrRd", rev = TRUE)
brks<-c(-0.01,seq(0.1,0.9,.1),1.01)

image.plot(aln_mat,axes=FALSE,col=cs,breaks=brks)
axis(1,at=(0:12)/12,usc,las=2)
axis(2,at=(0:13)/13,0:13)
box()
dev.off()
```
The map looks great: [GreenStr2Melanic.pdf](https://github.com/zgompert/TimemaFusion/files/7412533/GreenStr2Melanic.pdf)

Here is a table relating the chromosomes (16,151 was split between 9 and 13 and most of NA=0 was 14,101, which is the new 13)

| Melanic LG | Green Stripe chrom. |
|-----------:|--------------------:|
| 1 | 8483 |
| 2 | 14640 |
| 3 | 42935 |
| 4 | 42912 | 
| 5 | 18722 |
| 6 | 9928 |
| 7 | 10660 |
| 8 | 7748 |
| 9 | 16151 |
| 10 | 14160 |
| 11 | 12033 |
| 12 | 12380 |
| 0 | 14101 |

Old 0 is now and hereafter 13 (i.e., 14101 = LG 13) for new numbering.

* *T. knulli* genome versus the green striped *T. cristinae* genome. As with melanic versus green striped *T. cristinae*, I am using blast to detect homologous chromosomes between these species. 

I am again using the green striped blastable data base. I am working in `/uufs/chpc.utah.edu/common/home/u6000989/data/timema/hic_genomes/t_knulli`. There are 12 large scaffolds for *T. knulli* (not 13) consistent with karyotype data I think. I extracted these with `samtools` (version 1.12).

```{bash}
samtools faidx mod_hic_output.fasta -r knulliGenomeLGscafs.txt -o tknulli_chroms_hic_output.fasta
```
Then, I used `blastn` (version 2.11.0) for the local blast.

```{bash}
module load blast/2.11.0
cd /uufs/chpc.utah.edu/common/home/u6000989/data/timema/hic_genomes/t_knulli
blastn -db GreenGenome -evalue 1e-50 -perc_identity 92 -query tknulli_chroms_hic_output.fasta -outfmt 6 -num_threads 72 > Green2Knulli.txt
## column headers, see blastn -help for details
## qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
```

This did not seem to work so well. I am going to use a different approach.

* Trying *T. knulli* versus *T. cristinae* green striped with `cactus` as well. I am doing this because it is not clear (at this point) whether the approach using `blastn` will work with these more divergent genomes.

First, I am masking repeats with `RepateMasker` (version 4.0.7).

```{bash}
module load repeatmasker

#version 4.0.7
cd /uufs/chpc.utah.edu/common/home/u6000989/data/timema/hic_genomes/repeat_mask

## run repeat masker on each genome sequence, uses library from the 2020 Science paper 
## developed by Victor

RepeatMasker -s -e ncbi -xsmall -pa 24 -lib RepeatLibMergeCentroidsRM.lib /uufs/chpc.utah.edu/common/home/u6000989/data/timema/hic_genomes/t_crist/timema_cristinae_LGs_12Jun2019_lu3Hs.fasta
RepeatMasker -s -e ncbi -xsmall -pa 24 -lib RepeatLibMergeCentroidsRM.lib /uufs/chpc.utah.edu/common/home/u6000989/data/timema/hic_genomes/t_knulli/tknulli_chroms_hic_output.fasta
```

Next, I used `cactus` (version 1.0.0) to align the (repeat masked) *T. knulli* and green striped *T. cristinae* genomes. 

```{bash}
cd /scratch/general/lustre/cactusNp

module load cactus

cactus jobStore /uufs/chpc.utah.edu/common/home/u6000989/data/timema/hic_genomes/comp_aligns/cactusTimema_TcrGS_Tknul.txt cactusTcrGS_Tknul.hal  --maxCores 80 
```

Then, I used several HAL (Hierarchical Alignment) tools to summaize the alignments and extract synteny blocks. See [HAL](https://github.com/ComparativeGenomicsToolkit/hal) for a general description of these tools and [Krasheninnikova et al. 2020](https://academic.oup.com/gigascience/article/9/6/giaa047/5848161) for details about the `HalSynteny` algorithm. 

```{bash}
#!/bin/sh
## commands run thus far to summarize the hal alignment
module load cactus
## see https://github.com/ComparativeGenomicsToolkit/hal

~/source/hal/bin/halValidate cactusTcrGS_Tknul.hal
#File valid

~/source/hal/bin/halStats cactusTcrGS_Tknul.hal
#hal v2.1
#(t_cris_gs:0.01,t_knulli:0.01)Anc0;

#GenomeName, NumChildren, Length, NumSequences, NumTopSegments, NumBottomSegments
#Anc0, 2, 683142857, 42183, 0, 23649957
#t_cris_gs, 0, 1119939319, 13, 32708297, 0
#t_knulli, 0, 1145508109, 12, 36481639, 0

~/source/hal/bin/halSynteny --queryGenome t_knulli --targetGenome t_cris_gs cactusTcrGS_Tknul.hal out_synteny_knulli.psl
```

This generates a syntency file in psl format, for a description see [PSL](https://genome.ucsc.edu/FAQ/FAQformat.html#format2). The key columns am using are 1 = matches (number of matching bases), 10 = query name and 14 = target name. I summarized the matches between scaffolds/chromosomes to identify homologous chromosomes between *T. cristinae* and *T. knulli*. The results were unambiguous, see the table below and [SynPlotsKnulli.R](SynPlotsKnulli.R). Note, 1 and 3 were fused in *T. knulli* (both = scaffold 29).  [SynTcrTknul.pdf](https://github.com/zgompert/TimemaFusion/files/7500038/SynTcrTknul.pdf)

|T. christinae chrom. | Green Stripe scaf. | T. knulli scaf. |
|-----------:|--------------------:|-----------------------:|
| 1 | 8483 | 29 |
| 2 | 14640 | 813 |
| 3 | 42935 | 29 |
| 4 | 42912 | 6886 |
| 5 | 18722 | 6895 |
| 6 | 9928 | 6839 |
| 7 | 10660 | 934 |
| 8 | 7748 | 6852 |
| 9 | 16151 | 1305 |
| 10 | 14160 | 30 |
| 11 | 12033 | 500 |
| 12 | 12380 | 6840 |
| 13 | 14101 | 775 |

Next, I looked at patterns of colinearity for each pair of homolgous chromosomes (including 11). Structural variation between *T. cristinae* and *T. knulli*, see [SynPlotsKnulli.R](SynPlotsKnulli.R). Note that for this +- alignments are flipped. Structural variants are quite common and an inversion on Chromosome 11 ligns up with the PCA SV signal [AlnPlotsKnulTcr.pdf](https://github.com/zgompert/TimemaFusion/files/7671914/AlnPlotsKnulTcr.pdf).

I am trying an additional alignment with `mummer` (version 4.0.0rc1) to make sure I understand what the "strand" means. The program is described [here](https://github.com/mummer4/mummer). For this, I extracted just chromosome 11 from *T. knulli* (scaffold 500, no repeat masking) with `samtools` (version 1.12), and then ran the `nucmer` command from `mummer`.

```{bash}
samtools faidx ../t_knulli/tknulli_chroms_hic_output.fasta -o tknulli_chroms_hic_chrom11.fasta -r regionsq

## just chrom 11 from knulli against all lg chroms of cristinae
~/source/mummer-4.0.0rc1/nucmer -t 48 -p mumaln_tknulli_tcr tknulli_chroms_hic_chrom11.fasta /uufs/chpc.utah.edu/common/home/gompert-group1/data/timema/hic_genomes/t_crist/timema_cristinae_LGs_12Jun2019_lu3Hs.fasta.masked

~/source/mummer-4.0.0rc1/show-coords mumaln_tknulli_tcr.delta mumaln_tknulli_tcr.coords
```

As an additional comparison to understand structural variation in *T. knulli*, I am using `cactus` (version 1.0.0) to align our Hi-C *T. chumash* genome to the *T. knulli* genome. I am using a previously masked (with RepateMasker) version of the *T. chumash* genome.

```{bash}
module load cactus

## perform the alignment
cactus jobStore /uufs/chpc.utah.edu/common/home/u6000989/data/timema/hic_genomes/comp_aligns/cactusTimema_Tknul_Tchum.txt cactusTknul_Tchum.hal --maxCores 80 

## find synteny blocks
~/source/hal/bin/halSynteny --queryGenome t_chumash --targetGenome t_knulli cactusTknul_Tchum.hal out_synteny_KnulChum.psl
```
Summarized synteny blocks with [SynPlotsChumKnul.R](SynPlotsChumKnul.R). *T. chumash* has 10 chromsomes (this individual at least, I think Tanja's work suggest 11-13), with chromosome 1 representing four fused *T. cristinae* chrosomosomes (including 1 and 3 which were fused in *T. knulli*). Similar to the comparison with *T. cristinae*, structural variants are quite common and an inversion on Chromosome 11 ligns up with the PCA SV signal [AlnPlotsChumKnul.pdf](https://github.com/zgompert/TimemaFusion/files/7671919/AlnPlotsChumKnul.pdf).


|T. christinae chrom. | Green Stripe scaf. | T. knulli scaf. | T. chumash scaf. |
|-----------:|--------------------:|-----------------------:|-----------------------:|
| 1 | 8483 | 29 | 43 |
| 2 | 14640 | 813 | 1392 |
| 3 | 42935 | 29 | 43 |
| 4 | 42912 | 6886 | 43 |
| 5 | 18722 | 6895 | 56 |
| 6 | 9928 | 6839 | 1469 |
| 7 | 10660 | 934 | 1510 |
| 8 | 7748 | 6852 | 113 |
| 9 | 16151 | 1305 | 43 | 
| 10 | 14160 | 30 | 1213 |
| 11 | 12033 | 500 | 48 |
| 12 | 12380 | 6840 | 1403 |
| 13 | 14101 | 775 | 1308 |


Finally, mostly for completeness, I am using `cactus` (version 1.0.0) to align our Hi-C *T. chumash* genome to the green striped *T. cristinae* genome. I am using a previously masked (with RepateMasker) version of the *T. chumash* genome, which now just has the 10 scaffolds corresponding to chromsomes.

```{bash}

module load cactus

## whole genome alignment
cactus jobStore /uufs/chpc.utah.edu/common/home/u6000989/data/timema/hic_genomes/comp_aligns/cactusTimema_TcrGS_Tchum.txt cactusTcrGS_Tchum.hal  --maxCores 80 

## extract syneny blocks
~/source/hal/bin/halSynteny --queryGenome t_chum --targetGenome t_cris_gs cactusTcrGS_Tchum.hal out_synteny_CrisChum.psl
```

Summarized synteny blocks with [SynPlotsChumTcr.R](SynPlotsChumTcr.R). The synteny analysis is consistent with expectations from the alignments of both *T. chumash* and *T. cristinae* with *T. knulli* as summarized in the tables above. As expected, there is no evidence of a large inversion on chromosome 11 between *T. chumash* and *T. cristinae*, which further bolsters the evidence that we have an inversion in *T. knulli* relative to both of these species. See [AlnPlotsChumCris.pdf](https://github.com/zgompert/TimemaFusion/files/7804471/AlnPlotsChumCris.pdf).

Finally, finally (maybe), we now have a *T. podura* genome, so I am throwing that in the mix too. Here again, I first maksed repeats with `RepeatMasker` (version 4.0.7) and then performed the alignment (first to *T. knulli*) with `cactus` (version 1.0.0). I have a copy of the genome here: /uufs/chpc.utah.edu/common/home/gompert-group1/data/timema/hic_genomes/t_podura/.

```{bash}
#!/bin/sh

module load repeatmasker

#version 4.0.7
cd /uufs/chpc.utah.edu/common/home/u6000989/data/timema/hic_genomes/repeat_mask

RepeatMasker -s -e ncbi -xsmall -pa 24 -lib RepeatLibMergeCentroidsRM.lib /uufs/chpc.utah.edu/common/home/u6000989/data/timema/hic_genomes/t_podura/jasmine-cen2309-mb-hirise-gay8i__12-30-2021__hic_output.fasta

cd /scratch/general/lustre/cactusNp

module load cactus

cactus jobStore /uufs/chpc.utah.edu/common/home/u6000989/data/timema/hic_genomes/comp_aligns/cactusTimema_Tknul_Tpod.txt cactusTknul_Tpod.hal --maxCores 80 
```


## Sex chromosome

I used depth of coverage for *T. knulli* to identify the X sex chromosome (2 copies in females 1 in males). Presumably this is the same for *T. cristinae* but will check at some point. The depth data (from the redwood feeding experiment) is in `/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_fusion/genotypes_rw` (see `depthKnulli.txt`). The analysis is in [findSexChrom.R](findSexChrom.R). The X is chromosome 13 (as defined above), which mostly comprises parts of the genome that were not assigned to a linkage group (NA) in the old (pre Hi-C) melanic *T. cristinae* genome. See the sex-coverage plot.[SexChrom.pdf](https://github.com/zgompert/TimemaFusion/files/7541054/SexChrom.pdf)

## Nanopore data set

Linyi generated these data and called SVs. This data set comprises Oxford Nanopore MinION sequence data from 6 green *T. cristinae* from Refugio.

* Base calling with `guppy_basecaller`, which is part of the Guppy basecalling suite (version 4.2.2+effbaf8)
```{bash}
/uufs/chpc.utah.edu/common/home/u6033116/ont-guppy/bin/guppy_basecaller --input_path /uufs/chpc.utah.edu/common/home/gompert-group3/data/nanoporereads/TimemaL1/fast5/ --save_path /uufs/chpc.utah.edu/common/home/gompert-group3/data/nanoporereads/TimemaL1/HA_fastq --flowcell FLO-MIN106 --kit SQK-LSK109 -x auto --qscore_filtering --min_qscore 8
```

* De-multiplexing with `guppy_barcoder`, which is part of the Guppy basecalling suite (version 4.2.2+effbaf8)
```{bash}
/uufs/chpc.utah.edu/common/home/u6033116/ont-guppy/bin/guppy_barcoder --input_path /uufs/chpc.utah.edu/common/home/gompert-group3/data/nanoporereads/TimemaL1/HA_fastq/pass/ --save_path /uufs/chpc.utah.edu/common/home/gompert-group3/data/nanoporereads/TimemaL1/demux_fastq/ --config configuration.cfg --barcode_kits "EXP-NBD104" --trim_barcodes
```

* Sequences were aligned to the HiC green *T. cristinae* genome, `/uufs/chpc.utah.edu/common/home/gompert-group1/data/timema/hic_genomes/t_crist`. Need to add code from Linyi here still.

* Covert sam alignments to bam
```{bash}
for file in *.sam
do
    echo $file 
samtools view -S -b $file > ${file/.sam/.bam
```
* Call SVs with sniffles
```{bash}
for fi in *.bam; do
filename=$(awk 'FNR ==1{ print FILENAME}' "$fi"| awk -Fm '{print $1}')
sniffles -m $fi -v "$filename"gt.vcf --Ivcf TM_merged_SURVIVOR_1kbpdist_typesavenew.vcf
done

}
```

* Translocation file from Linyi 

[Timema_TRA.txt](https://github.com/zgompert/TimemaFusion/files/7384224/Timema_TRA.txt)

* Code and ideas for visualizing translocations

[SV annotation](https://bioconductor.org/packages/devel/bioc/vignettes/StructuralVariantAnnotation/inst/doc/vignettes.html)



* **NOTE** Might want to try *de novo* assembly of nanopore reads instead, see [canu](https://github.com/marbl/canu). Also, evidence from TRA transition calls doesn't really support merger of LG 8 and 11 (7748 and 12033), or of 42912 and 42934, which show highest interchromosomal LD in Refugio. TRA call might not be that useful. Instead use these data for *de novo* assembly and maybe for identifying other SVs (inversions and deletions) relative to Hwy154 (nice because both based on green genome).

## Host-associated genetic differentiation for 2017 Nature EE populations

To motivate the analysis of RW adaptation in *T. knulli*, I want to quantify host-plant associated differentiation across other northern clade *Timema* species. For this, I am using the GBS data from [Chaturvedi et al. 2021-preprint](https://www.researchsquare.com/article/rs-923547/v1), which were re-analyzed by Sam for the climate-adaptation paper [Riesch et al. 2017](https://www.nature.com/articles/s41559-017-0082). I will get processing details from that paper. I am starting from her genotype likelihood files. These, along with a file with IDs are in `/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_fusion/host_diff`. This includes data for *T. californicum*, *T. knulli*, *T. landelsensis* and *T. poppensis*. 

I have split the data by population and host with [splitPops.pl](splitPops.pl).

Here are the sample sizes, sites and hosts for the populations I will use:

| Species |Pop. 1 | Pop. 2 | N1 | N2 | No. SNPs |
|--------:|------:|-------:|---:|---:|---------:|
| *T. californicum* | SM on M | SM on Q| 17 | 20 | 7858|
| *T. knulli* | BCE on RW | BCWP on C| 15 | 12 | 1139|
| *T. knulli* | BCTUR on C | BCTUR on P| 17 | 16 | 1139|
| *T. landelsensis* | BCBOG on C | BCBOG on Q| 23 | 20 | 8548|
| *T. landelsensis* | BCSUM on C | BCSUM on Q| 20 | 11 | 8548|
| *T. poppensis* | TBARN on DF | TBARN on RW| 20 | 20 | 7157|

I calculated Fst for each SNP and summarized Fst at the linkage group level (for the old *T. cristinae* LGs from the melanic genome). See [hostFst.R](hostFst.R). The signal for *T. knulli* RW vs C is clear (elevated Fst on LG 11). Nothing else really stands out, except a weaker signal on LG 8 likely associated with color in *T. californicum*.

## Alignment, variant calling and filtering for GBS data

* **rw_plus** data set; includes *T. knulli* and *T. petita*.

1. DNA sequence alignment with `bwa`; used the the *T. knulli* genome.

2. Compress, sort and index alignments with `samtools`.

3. Variant calling with `samtools` (version 1.5) and `bcftools` (version 1.6).

For *T. knulli*:

```{bash}
module load samtools/1.5
module load bcftools
## samtools 1.5
## bcftools 1.6
cd /uufs/chpc.utah.edu/common/home/gompert-group3/data/timema_clines_rw_SV/align_rw_plus
samtools mpileup -b bams_knulli -C 50 -d 500 -f /uufs/chpc.utah.edu/common/home/u6000989/data/timema/hic_genomes/t_knulli/mod_hic_output.fasta -q 20 -Q 30 -I -g -u -t DP,AD,ADF,ADR -o tcr_rw_knulli_variants.bcf
bcftools call -v -c -p 0.01 -O v -o tcr_rw_knulli_variants.vcf tcr_rw_knulli_variants.bcf
```

For *T. petita*:
```{bash}
module load samtools/1.5
module load bcftools
## samtools 1.5
## bcftools 1.6
cd /uufs/chpc.utah.edu/common/home/gompert-group3/data/timema_clines_rw_SV/align_rw_plus
samtools mpileup -b bams_petita -C 50 -d 500 -f /uufs/chpc.utah.edu/common/home/u6000989/data/timema/hic_genomes/t_knulli/mod_hic_output.fasta -q 20 -Q 30 -I -g -u -t DP,AD,ADF,ADR -o tcr_rw_petita_variants.bcf
bcftools call -v -c -p 0.01 -O v -o tcr_rw_petita_variants.vcf tcr_rw_petita_variants.bcf
```

For *T. knulli* and *T. petita* combined (this is to have a common set of variants for a phylogenetic analysis to date the perform SV locus)

```{bash}
module load samtools/1.5
module load bcftools
## samtools 1.5
## bcftools 1.6
cd /uufs/chpc.utah.edu/common/home/gompert-group3/data/timema_clines_rw_SV/align_rw_plus
samtools mpileup -b bams_comb -C 50 -d 500 -f /uufs/chpc.utah.edu/common/home/u6000989/data/timema/hic_genomes/t_knulli/mod_hic_output.fasta -q 20 -Q 30 -I -g -u -t DP,AD,ADF,ADR -o tcr_rw_comb_variants.bcf
bcftools call -v -c -p 0.01 -O v -o tcr_rw_comb_variants.vcf tcr_rw_comb_variants.bcf
```

* Variant filtering with `vcfFilterKnullipl`, `vcfFilterPetita.pl` and `filterSomeMore.pl`.

I filtered based on the same criteria for both species/data sets: 2X minimum coverage per individual, a minimum of 10 reads supporting the alternative allele, Mann-Whittney P values for BQ, MQ and read position rank-sum tests > 0.005, a minimum ratio of variant confidence to non-reference read depth of 2, a minimum mapping quality of 30, no more than 20% of individuals with missing data, only bi-allelic SNPs, and coverage not > 3 SDs of the mean coverage (at the SNP level).

This left me with **64,650** SNPs for *T. knulli* (N = 138 individuals) and **32,859** SNPs for *T. petita* (N = 69 individuals).

* Lastly, vcf files were converted to genotype likelihood format.

```{bash}
perl vcf2glSamt.pl 0.0 morefilter_filtered2x_tcr_rw_knulli_variants.vcf
#Number of loci: 64650; number of individuals 138
perl vcf2glSamt.pl 0.0 morefilter_filtered2x_tcr_rw_petita_variants.vcf 
#Number of loci: 32859; number of individuals 69
```

## Alignment, variant calling and filtering for mate-pair data

## Genotype inference for *T. knulli* and *T. petita* (redwood project)

Working in `/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_fusion/genotypes_rw` with gl files linked from `/uufs/chpc.utah.edu/common/home/gompert-group3/data/timema_clines_rw_SV/variants_rw_plus`.

I estimated genotypes using `entropy` (version 1.2). Estimates were based on k = 2,3 and used starting values from LDA of PC scores.

* Generate initial values for admixture proportions.

Calculate simple genotype point estimates.
```{bash}
perl gl2genest.pl af_filtered_tknulli_variants.txt filtered_tknulli_variants.gl
perl gl2genest.pl af_filtered_tpetita_variants.txt filtered_tpetita_variants.gl
```
Generate initial values for admixture proportions.
```{R}
## read genotype estimates, i.e., g = 0 * Pr(g=0) + 1 * Pr(g=1) + 2 * Pr(g=2)
## row = locus, column = individual
L<-64650
N<-138
g<-matrix(scan("pntest_filtered_tknulli_variants.txt",n=N*L,sep=" "),nrow=L,ncol=N,byrow=T)

## pca on the genotype covariance matrix
pcgcov<-prcomp(x=t(g),center=TRUE,scale=FALSE)

## kmeans and lda
library(MASS)
k2<-kmeans(pcgcov$x[,1:5],2,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k3<-kmeans(pcgcov$x[,1:5],3,iter.max=10,nstart=10,algorithm="Hartigan-Wong")

ldak2<-lda(x=pcgcov$x[,1:5],grouping=k2$cluster,CV=TRUE)
ldak3<-lda(x=pcgcov$x[,1:5],grouping=k3$cluster,CV=TRUE)

write.table(round(ldak2$posterior,5),file="ldak2_knulli.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak3$posterior,5),file="ldak3_knulli.txt",quote=F,row.names=F,col.names=F)

save(list=ls(),file="init_knulli.rdat")

L<-32859
N<-69
g<-matrix(scan("pntest_filtered_tpetita_variants.txt",n=N*L,sep=" "),nrow=L,ncol=N,byrow=T)

## pca on the genotype covariance matrix
pcgcov<-prcomp(x=t(g),center=TRUE,scale=FALSE)

## kmeans and lda
library(MASS)
k2<-kmeans(pcgcov$x[,1:5],2,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k3<-kmeans(pcgcov$x[,1:5],3,iter.max=10,nstart=10,algorithm="Hartigan-Wong")

ldak2<-lda(x=pcgcov$x[,1:5],grouping=k2$cluster,CV=TRUE)
ldak3<-lda(x=pcgcov$x[,1:5],grouping=k3$cluster,CV=TRUE)

write.table(round(ldak2$posterior,5),file="ldak2_petita.itxt",quote=F,row.names=F,col.names=F)
write.table(round(ldak3$posterior,5),file="ldak3_petita.txt",quote=F,row.names=F,col.names=F)

save(list=ls(),file="init_petita.rdat")

## when you run entropy use provide the input values as, e.g., -q ldak2.txt
## also set -s to something like 50
```
Bayesian genotype estimates were then obtained running entropy as follows (input files = `filtered_tknulli_variants.gl` and `filtered_tpetita_variants.gl`):

```{bash}
#!/bin/sh 
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --account=gompert
#SBATCH --partition=notchpeak
#SBATCH --job-name=entropy
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load gsl
module load hdf5

cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_fusion/genotypes_rw

perl RunEntropyFork.pl filtered_*gl
```

```{perl}
#!/usr/bin/perl
#
# run entropy jobs
#

use Parallel::ForkManager;
my $max = 20;
my $pm = Parallel::ForkManager->new($max);

my $odir = "/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_fusion/genotypes_rw/Entropy/";
my $base = "/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_fusion/genotypes_rw/";

foreach $in (@ARGV){
	$in =~ m/^([a-zA-Z0-9_\-]+)/;
	$dat = $1;
	$in =~ m/_t([a-z]+)_/;
	$sp = $1;
	sleep(2);
	foreach $k (2..3){
		foreach $ch (0..4){
		$pm->start and next; ## fork
		system "entropy -i $in -l 8000 -b 5000 -t 3 -k $k -Q 0 -s 50 -q $base"."ldak$k"."_$sp.txt -o $odir"."out_$dat"."_k$k"."_ch$ch".".hdf5 -w 0 -m 1\n";
		$pm->finish;
		}
	}
}

$pm->wait_all_children;
```
Last, I used `estpost.entropy` (version 2.0) to summarize the posteriores and specifically obtain point estimates of genotypes averaging over the 5 chains and 2 values of *k*.

```{bash}
estpost.entropy -p gprob -s 0 -w 0 *hdf5 -o G_tcr_refugio.txt
#parameter dimensions for gprob: loci = 64650, ind = 138, genotypes = 3, chains = 10
estpost.entropy -p gprob -s 0 -w 0 out_filtered_tpetita_variants_k*hdf5 -o G_tpetita.txt
#parameter dimensions for gprob: loci = 32859, ind = 69, genotypes = 3, chains = 10
```

## Multilocus mapping from *T. knulli* and *T. petita* experiments

Old data and analyses are in `/uufs/chpc.utah.edu/common/home/u6000989/projects/timema_confiers/redwood_gwa/`. My original analyses are described on the legacy lab site [Timema redwood](https://sites.google.com/site/gompertlabnotes/home/researcher-pages/zach-gompert/timema/timema-redwood). 

The goal is to map performance (weight, weight change between 15 and 21 days and survival) on *Ceanothus* and redwoods.

* First, I formatted the genetic and phenotypic data from the experiment for mapping with gemma. This includes splitting the data set by treatment (host) and removing the effects of initial weight/stage on the performance traits. I am using these five traits:

1. 15 day weight, control sex and stage
2. 21 day weight, control sex and stage
3. Survival
4. 21-15 day weight change, control sex and stage
5. 21-15 day weight change, control sex

See [formatPhenoGeno.R](formatPhenoGeno.R).

* Add SNP names and placeholder alleles IDs (don't really matter) to the genotype files for `gemma`.

```{bash}
perl mkGenoFile.pl geno_*
```

```{perl}
#!/usr/bin/perl
#
# formats the geno file
#

## get actual SNP ids
open(IN, "Snps.txt") or die "failed to read SNPs file\n";
while(<IN>){
        chomp;
        push(@snps,$_);
}
close(IN);

foreach $in (@ARGV){ ## geno files

        ## read and write geno file
        open(IN, $in) or die "failed to read the genotype file\n";
        $out = "mod_$in";
        open(OUT, "> $out") or die "failed to write\n";
        $i = 0;
        while(<IN>){
                chomp;
                s/ /, /g or die "failed space sub\n";
                print OUT "$snps[$i], A, T, $_\n";
                $i++;
        }
        close(IN);
        close(OUT);
}
```

* Run the BSLMM in `gemma` (Version 0.95alpha), 10 chains per phenotype, both treatments, full and subset data sets without BCTURN and all five performance traits.

```{bash}
module load gemma
# Version 0.95alpha

cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_fusion/redwood_gwa

perl RunGemmaFork.pl mod_geno*
```

```{perl}
#!/usr/bin/perl
#
# run gemma jobs
#

use Parallel::ForkManager;
my $max = 48;
my $pm = Parallel::ForkManager->new($max);

foreach $g (@ARGV){
	$ph = $g;
	$ph =~ s/mod_g/ph/;
	$base = $g;
	$base =~ s/mod_geno_//;
	$base =~ s/\.txt//;
	foreach $ch (0..9){
		system "sleep 2\n";
		foreach $trait (1..5){
			$pm->start and next; ## fork
			$out = "o_tknulli_$base"."_ph$trait"."_ch$ch";
			system "gemma -g $g -p $ph -bslmm 1 -n $trait -maf 0.0 -w 200000 -s 1000000 -o $out\n";
			$pm->finish;
		}
	}
}

$pm->wait_all_children;
```

* Next, I summarized the relevant posteriors for mapping.

Posterior estimates of hyperparameters (see [calpost.pl](calpost.pl)).
```{bash}
perl calpost.pl o_tknulli_RW_ph*ch0.hyp.txt
perl calpost.pl o_tknulli_RW_sub_ph*ch0.hyp.txt
perl calpost.pl o_tknulli_C_ph*ch0.hyp.txt
perl calpost.pl o_tknulli_C_sub_ph*ch0.hyp.txt
```

Mostly these seem reasonable, though survival (phenotype 3) seems too high (probably not reliable).

Extracted posterior inclusion probabilities and made some plots of PIPs, no real evidence that scaffol 500 (the SV scaffold) pops out (see [grabPips.pl](grabPips.pl)).
```{bash}
perl grabPips.pl o_tknulli_*ch0.param.txt
```
```{R}
## summarize PIPs for knulli
pf<-list.files(pattern="pip") ## 20 files, C, Csub, RW, RW sub, 5 pheno each
L<-64650
K<-20

pips<-matrix(NA,nrow=L,ncol=K)
for(i in 1:K){
        pips[,i]<-scan(pf[i])
}

## on C
par(mfrow=c(2,5))
for(i in 1:10){
        plot(pips[,i])
}

## on RW
par(mfrow=c(2,5))
for(i in 11:20){
        plot(pips[,i])
}

snps<-read.table("../output/snpList.txt",header=FALSE)

sc500<-snps[,1]==500


cs<-c("darkgray","red")
cl<-1.5;cm<-1.5
titls<-rep(c("w15d","w21d","surv","dW","dW2"),4)

pdf("KnulliPipsC.pdf",width=10,height=12)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,1))
for(i in 1:5){
        a<-which(pips[,i] > 0.001) 
        plot((1:L)[a],pips[a,i],pch=19,col=cs[sc500+1][a],xlab="SNP number",ylab="PIP",cex.lab=cl)
        title(main=titls[i],cex.main=cm)
}
dev.off()
## same for other combinations, see summarizePips.R
```
I also tried calling "genotypes" for the "SV" on scaffold 500 (PCA and k-means clustering). I then asked whether SV genotype is associated with performnance (see [lmSV.R](lmSV.R)). There is some evidence of a negative association with 15 d weight on RW, 15 and 21 d weight on C, and a positive association on C. Interesting but need to think more.

```{R}
## R script to test SV/fusion genotype effect on performance
load("gp.rdat")

## scaffold 500
a<-which(snps[,1]==500)

## pca RW without BCTURN
pc<-prcomp(t(g_RW_sub[a,]),center=TRUE,scale=FALSE)
ko<-kmeans(pc$x[,1],centers=3)
gen<-ko$cluster ## these are sorted

## lm fits, RW only 15 d weight significant
summary(lm(gemma_phRW_sub[,1] ~ gen))

```

## Gene flow-selection balance for the SV locus in *T. knulli*

See `/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_fusion/TknulliGF`.

Estimating gene flow for BCE C, BCE RW and BCTURN (allopatric on C) using a Bayesian F-model, [beta_ismod.stan](beta_ismod.stan) and [MigStan.R](MigStan.R). Simulations showing directional selection gene flow balance [sim_sm_balance.R](sim_sm_balance.R). My plan is to develop an ABC model that takes advantage of Nm for alloptaric vs sympatric C x RW.

## LD based Ne estiamtes to parameteriz ABC model

Method to estimate delta follows [Weir 1979](https://www.jstor.org/stable/2529947?seq=1#metadata_info_tab_contents) and [Zaykin 2004](https://onlinelibrary.wiley.com/doi/epdf/10.1002/gepi.20015). From there, I wrote my own R script for LD-based estimationg of Ne following [Waples and Do 2008](https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1755-0998.2007.02061.x). A distinct angle of my approach was to subsample independent pairs of SNP loci on different chromosomes. The main functions for the LD-based inference of Ne are in [ldNe.R](ldNe.R) and the analysis is and results are in [knulliEstLdNe.R](knulliEstLdNe.R). Posterior means for Ne (not 2Ne) were 23.2 for BCTURN C, 125.0 for BCE RW and 75.7 for BCE C.

## ABC model of D/B selection gene flow balance

Simulations for fitting a gene flow-selection model were run the a program written in C++. See [main.C](main.C), [func.C](func.C) and [kbsel.H](kbsel.H).

Ran simulations with 3 sets of priors, based either on the contemporary estimates of Ne from LD, or wider priors that I think better reflect our actual undertainty. 25,000,000 simulations of each.

```{perl}
#!/usr/bin/perl
#
# fork run of ABC modelx 
#


use Parallel::ForkManager;
my $max = 50;
my $pm = Parallel::ForkManager->new($max);

foreach $k (0..99){
        sleep(2);
        $pm->start and next; ## fork
        system "~/bin/kbsel -a 50 -b 50 -c 50 -A 1000 -B 1000 -C 1000 -n 250000 -o out_wide_$k".".txt\n";
        system "~/bin/kbsel -a 100 -b 100 -c 100 -A 2000 -B 2000 -C 2000 -n 250000 -o out_hi_$k".".txt\n";
        system "~/bin/kbsel -a 22 -b 96 -c 72 -A 24 -B 137 -C 79 -n 250000 -o out_est_$k".".txt\n";
        $pm->finish;
}

$pm->wait_all_children;
```

I then used the `abc` (2.1) package in `R` (version 4.0.2) to summarize posteriors. I used simple rejection for model choice. The model with balancing selection on both hosts was best under all priors. Model-averaging doesn't make sense across models here (different parameters), so I inferred selection based on the balancing selection model. I used the log transformation of the selection coefficients and ridge regression for parameter adjustment. See [SummarizeABC.R](SummarizeABC.R).


## Trying to determine the nature of the SV locus in *T. knulli*

The patterns of genetic variation on chromosome 11 (scaffold 500) in *T. knulli* indicate a large region of reduced recombination, which is most likely some form of structural variant. I looked at patterns of coverage (for SNPs) and LD along that chromosome and by SV genotype (i.e., PC1 cluster) within *T. knulli* with hopes of shedding light on the nature of the SV. There were some weak signals associated with coverage (maybe a indel is involved near one edge) with more compelling patterns in terms of LD (see [LDdifplot.pdf](https://github.com/zgompert/TimemaFusion/files/7507542/LDdifplot.pdf)
[LDplot.pdf](https://github.com/zgompert/TimemaFusion/files/7507543/LDplot.pdf)). In brief, LD is elevated in the C allele relative to the RW allele across the SV region, and especially at the boundaries. This is consistent with selection favoring the C allele or with less recombinatino within C than within RW, but still doesn't demonstrate conclusively the nature of the SV. See [knulliDepthLD.R](knulliDepthLD.R).

I used the eigenvalues from a PCA in 100 SNP windows along chromosome 11 (scaffold 500) as an alternative approach to determine the SV boundaries (but not types). The eigenvalues increase when you hit the SV as the first PC explains more of the total variation (this excludes BCTURN to avoid confounding from actual structure). I then fit a HMM to precisely define the boundaries. This with done in `R` with `HiddenMarkov` version (1.8.13). See [defineBounds.R](defineBounds.R) for details. The HMM identified a single, continuous elevated eigenvalue region from 13,093,370 43,606,674. I will use this as the SV bounds for now at least.

With the latest alignments between *T. knulli* and (i) *T. cristiane* or *T. chumash*, I think we can be confident that this is an inversion. At minimum, the individual we sequenced is inversted relative to these species with the inversion corresponding precisely to the PCA SV signal. I still need to figure out whether the inverted allele is the RW or C allele and try to date the inversion.

## Divergence time dating for the SV locus alleles

My plan is to determine the divergence time for the two SV alleles (C vs. RW) for the perform locus. This can be done in a phylogeneitc context. We have a tree from Victor, described in [Riesch_et_al-2017](https://github.com/zgompert/TimemaFusion/files/7738281/Riesch_et_al-2017-Nature_Ecology_.26_Evolution.pdf). Doro used this tree for divergence time data of the *Mel-Stripe* [Lindtke et al. 2017](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14280). The key is to use the callibrations from Victor's tree to callibrate a tree with *T. knullia* SV alleles and *T. petita*. See the DRYAD code from [Doro](https://datadryad.org/stash/dataset/doi:10.5061/dryad.jt644) and [Victor](https://datadryad.org/stash/dataset/doi:10.5061/dryad.nq67q), along withe Doro's [supplemental material](https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fmec.14280&file=mec14280-sup-0001-Supinfo.pdf).


## LD for refugio versus hwy154

Working from `ld_refugio` and `ld_hwy154` within `/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_fusion`. 

* Used `entropy` (version 1.2) to estimate genotypes for both populations.

Perl submission script for Refugio:
```{perl}
#!/usr/bin/perl
#
# run entropy jobs
#

use Parallel::ForkManager;
my $max = 20;
my $pm = Parallel::ForkManager->new($max);

my $odir = "/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_fusion/ld_refugio/Entropy";
my $base = "/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_fusion/ld_refugio/";

foreach $in (@ARGV){
	$in =~ m/^([a-zA-Z0-9_\-]+)/;
	$dat = $1;
	foreach $k (2..3){
		foreach $ch (0..4){
		$pm->start and next; ## fork
		system "entropy -i $in -l 8000 -b 5000 -t 3 -k $k -Q 0 -s 50 -q $base"."ldak$k".".txt -o $odir"."out_$dat"."_k$k"."_ch$ch".".hdf5 -w 0 -m 1\n";
		$pm->finish;
		}
	}
}

$pm->wait_all_children;
```
And then conversion to genotype point estimateas (posteior over k = 2 and 3).

```{bash}
estpost.entropy -p gprob -s 0 -w 0 *hdf5 -o G_tcr_refugio.txt
```
Identical procedure for genotype estimation for Hwy154.

* LD between chromosomes (13 big scaffolds).

Next I calculated with and between chromosome (13 big scaffolds) LD for Hwy154 and Refugio. The key results and files are in `/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_fusion/ld_hwy154` and  `/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_fusion/ld_Refugio`. Here is the R script (InterChromLD.R).

```{R}
library(data.table)
library(RColorBrewer)
## first hwy154

## read genotype estimates from entropy
G<-as.matrix(fread("G_tcr_hwy154.txt",sep=",",header=FALSE))

dim(G)
#[1]   596 12616
# 596 inds. 12616 SNPs

## PCA mostly for fun, see how much structure
o<-prcomp(G,center=TRUE)
## PC 1 ~2%, clearly two clusters

## get snp data, scaffold and position
snps<-as.matrix(read.table("SNPs_scaf_pos.txt",header=FALSE))
scafs<-scan("bigScafs.txt")

## split G into list by LG=scaffold
G_LG<-vector("list",13)
for(k in 1:13){
	x<-which(snps[,1]==scafs[k])
	G_LG[[k]]<-G[,x]
}

## get mean correlation between LGs
Chrom_cors<-matrix(NA,nrow=(13*12)/2,ncol=3)
k<-1
for(a in 1:12){for(b in (a+1):13){
	Chrom_cors[k,1]<-scafs[a]
	Chrom_cors[k,2]<-scafs[b]
	o<-cor(G_LG[[a]],G_LG[[b]])
	Chrom_cors[k,3]<-mean(o^2)	
	k<-k+1
}}

## now Refugio

## read genotype estimates from entropy
G_r<-as.matrix(fread("../ld_refugio/G_tcr_refugio.txt",sep=",",header=FALSE))

dim(G_r)
#[1]   238 74658
# 238 inds. 74658 SNPs

## PCA mostly for fun, see how much structure
o<-prcomp(G_r,center=TRUE)
## PC 1 ~2%, gradient on 1, a few oddballs on 2

## get snp data, scaffold and position
snps_r<-as.matrix(read.table("../ld_refugio/SNPs_scaf_pos.txt",header=FALSE))

## split G into list by LG=scaffold
G_LG_r<-vector("list",13)
for(k in 1:13){
	x<-which(snps_r[,1]==scafs[k])
	G_LG_r[[k]]<-G_r[,x]
}

## get mean correlation between LGs
Chrom_cors_r<-matrix(NA,nrow=(13*12)/2,ncol=3)
k<-1
for(a in 1:12){for(b in (a+1):13){
	cat(paste(a,b,"\n"))
	Chrom_cors_r[k,1]<-scafs[a]
	Chrom_cors_r[k,2]<-scafs[b]
	o<-cor(G_LG_r[[a]],G_LG_r[[b]])
	Chrom_cors_r[k,3]<-mean(o^2)	
	k<-k+1
}}

## comparison
plot(Chrom_cors[,3],Chrom_cors_r[,3],pch=19,xlab="Mean LD (r2) Hwy154",ylab="Mean LD (r2) Refugio")
o<-lm(Chrom_cors_r[,3] ~ Chrom_cors[,3])
abline(o$coefficients)

## this one stands out, make a plot
Chrom_cors_r[which(Chrom_cors_r[,3] > 0.0047),]
## 42912 vs 42935
a<-9 ## 42912
b<-10 ## 42935
cc_9x10_h<-cor(G_LG[[a]],G_LG[[b]])
cc_9x10_r<-cor(G_LG_r[[a]],G_LG_r[[b]])
cc_9x9_h<-cor(G_LG[[a]],G_LG[[a]])
cc_10x10_h<-cor(G_LG[[b]],G_LG[[b]])
cc_9x9_r<-cor(G_LG_r[[a]],G_LG_r[[a]])
cc_10x10_r<-cor(G_LG_r[[b]],G_LG_r[[b]])

bnds<-c(-1,0.005,0.01,0.05,0.1,0.5,1)
bnds<-c(-1,0.1,0.5,1)
cs<-brewer.pal(n=6,"Reds")[c(1,5,6)]

image(cc_9x10_h^2,col=cs,breaks=bnds)
image(cc_9x9_h^2,col=cs,breaks=bnds)
image(cc_10x10_h^2,col=cs,breaks=bnds)
image(cc_9x10_r^2,col=cs,breaks=bnds)
image(cc_9x9_r^2,col=cs,breaks=bnds)
image(cc_10x10_r^2,col=cs,breaks=bnds)


for(a in 1:12){for(b in (a+1):13){
	nm<-paste("ld_",scafs[a],"_",scafs[b],".png",sep="")	
	oh_11<-cor(G_LG[[a]],G_LG[[a]])^2
	oh_22<-cor(G_LG[[b]],G_LG[[b]])^2
	oh_12<-cor(G_LG[[a]],G_LG[[b]])^2
	or_11<-cor(G_LG_r[[a]],G_LG_r[[a]])^2
	or_22<-cor(G_LG_r[[b]],G_LG_r[[b]])^2
	or_12<-cor(G_LG_r[[a]],G_LG_r[[b]])^2
	png(nm,width=12000,height=18000,units="px")
	par(mfrow=c(3,2))
	image(oh_11,col=cs,breaks=bnds)
	image(or_11,col=cs,breaks=bnds)
	image(oh_22,col=cs,breaks=bnds)
	image(or_22,col=cs,breaks=bnds)
	image(oh_12,col=cs,breaks=bnds)
	image(or_12,col=cs,breaks=bnds)
	dev.off()
}}
save(list=ls(),file="LD.rdat")
}}
```


Patterns of LD most consistent with a fusion of 42912 and 42935 for Refugio. Also the inversion on LG8=7748 is really apparent in the LD plots. At this point, there is nothing obvioulsy pointing to a fusion invovling 7748 (LG8), but let's see what the nanopore data show.
