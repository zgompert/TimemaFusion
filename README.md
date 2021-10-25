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
