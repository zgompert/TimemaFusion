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

From past alignments (details to come) Green 12033 is LG11 (which we are thinking about for the redwood stuff) and Green 7748 = LG8.

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

* Identifying pairs of chromosomes with translocations

The 13 (unambiguous) large scaffold (> 10 Mbps) from the green genome are listed, along with their sizes in bps, in `/uufs/chpc.utah.edu/common/home/gompert-group1/data/timema/hic_genomes/t_crist/greenGenomeLGscafs.txt`. I used this to extract these scaffolds from the genome with `samtools` (version 1.12)

```{bash}
samtools faidx timema_cristinae_12Jun2019_lu3Hs.fasta -r greenLGRegions.txt -o timema_cristinae_LGs_12Jun2019_lu3Hs.fasta
```
Next I created a blastable data base to match scaffolds from the melanic genome to the green genome. This was done with `blast` (version 2.11.0).

```{bash}
makeblastdb -in timema_cristinae_LGs_12Jun2019_lu3Hs.fasta -dbtype nucl -parse_seqids -out GreenGenome -title "Tcr Green genome chrom. scafs."
```

Trying to use `blastn` to match chromosomes

```{bash}
blastn -db GreenGenome -evalue 1e-50 -perc_identity 92 -query ../../tcrDovetail/version3/mod_map_timema_06Jun2016_RvNkF702.fasta -outfmt 6 -num_threads 48 > Green2Melanic.txt
```

* **NOTE** Might want to try *de novo* assembly of nanopore reads instead, see [canu](https://github.com/marbl/canu). Also, evidence from TRA transition calls doesn't really support merger of LG 8 and 11 (7748 and 12033), or of 42912 and 42934, which show highest interchromosomal LD in Refugio. TRA call might not be that useful. Instead use these data for *de novo* assembly and maybe for identifying other SVs (inversions and deletions) relative to Hwy154 (nice because both based on green genome).


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

* Variant filtering with `vcfFilterKnullipl`, `vcfFilterPetita.pl` and `filterSomeMore.pl`.

I filtered based on the same criteria for both species/data sets: 2X minimum coverage per individual, a minimum of 10 reads supporting the alternative allele, Mann-Whittney P values for BQ, MQ and read position rank-sum tests > 0.005, a minimum ratio of variant confidence to non-reference read depth of 2, a minimum mapping quality of 30, no more than 20% of individuals with missing data, only bi-allelic SNPs, and coverage not > 3 SDs of the mean coverage (at the SNP level).

This left me with **64,650** SNPs for *T. knulli* (N = 138 individuals) and **32,859** SNPs for *T. petita* (N = 69 individuals).


## Alignment, variant calling and filtering for mate-pair data

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
