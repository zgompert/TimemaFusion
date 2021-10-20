# TimemaFusion

Characterizing chromosomal fusion(s) in *Timema* stick insects

## Data sets

Most of the data sets for this project are in `/uufs/chpc.utah.edu/common/home/gompert-group3/data/timema_clines_rw_SV/`; the nanopore data for *T. cristinae* from Refugio are in `/uufs/chpc.utah.edu/common/home/gompert-group3/data/nanoporereads/TimemaL1`.

0. Whole genomes for compariative genome alignment
1. **hwy154** = GBS data
2. **fha_mapping_sample** = GBS data
3. **n1_doro** = GBS data
4. **refugio** and **refugio_2019** = GBS data
5. **rw_plus** = GBS data
6. **refugio_nanoper** = nanpore data from 6 *Timema cristinae* from Refugion, could provide direct evidence of a fusion on Refugio. 
7. **matepair_SVs** = matepair data form Kay's paper, could provide direct evidence of a fusion on Refugio.

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



## Alignment, variant calling and filtering for GBS data

## Alignment, variant calling and filtering for WGS data

## Alignment, variant calling and filtering for mate-pair data

