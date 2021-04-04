# AlleleSeq2

## Generate personal genomes, STAR indices and other helper files:



```
wgs_bam_file=genotyping_HTS/ENCFF200XCY.bam


make -f ~/bin/AlleleSeq2/makePersonalGenome.mk \
        N_THREADS=8 \
        ALIGNER=STAR \
        VCF_SAMPLE_ID=sge_Aug_encodev2_2_local \
        REFGENOME_VERSION=GRCh38 \
        OUTPUT_DIR=pgenome_ENC-003 \
        FILE_PATH_BAM=${wgs_bam_file} \
        FILE_PATH_VCF=enc003.spliced.scrubbed.vcf 
```



## (1) Call ASE hetSNVs from a single sample ###

```
r1=ENCFF337ZBN.fastq.gz
r2=ENCFF481IQE.fastq.gz  

pgenome=../../../../pgenomes_20191127/pgenome_ENC-003

make -f ~/bin/AlleleSeq2/PIPELINE.mk \
        PGENOME_DIR=${pgenome} \
        REFGENOME_VERSION=GRCh38 \
        NTHR=8 \
        READS_R1=${r1} \
        READS_R2=${r2} \
        STAR_readFIlesCommand=zcat \
        ALIGNMENT_MODE=ASE \
        RM_DUPLICATE_READS=on \
        FDR_CUTOFF=0.10 \
	VCF_SAMPLE_ID=sge_Aug_encodev2_2_local
```

#### the main output containing ASE hetSNVs is 
#### ENCFF337ZBN_ENCFF481IQE_interestingHets.FDR-0.10.betabinom.chrs1-22.6-tot_0-min_cnt.tsv





## (2) Pool two replicates, if available ###

```
pgenome=../../../../pgenomes_20191127/pgenome_ENC-003

make -f ~/bin/AlleleSeq2/PIPELINE_aggregated_counts.mk \
        PGENOME_DIR=${pgenome} \
        INPUT_UNIQ_READS_PILEUP_FILES="../ENCSR238ZZD_ENCFF719MSG_1_ENCFF120MML_2_1_1/ENCFF719MSG_ENCFF120MML_hap1_uniqreads.mpileup    ../ENCSR238ZZD_ENCFF719MSG_1_ENCFF120MML_2_1_1/ENCFF719MSG_ENCFF120MML_hap2_uniqreads.mpileup ../ENCSR238ZZD_ENCFF337ZBN_1_ENCFF481IQE_2_1_1/ENCFF337ZBN_ENCFF481IQE_hap1_uniqreads.mpileup ../ENCSR238ZZD_ENCFF337ZBN_1_ENCFF481IQE_2_1_1/ENCFF337ZBN_ENCFF481IQE_hap2_uniqreads.mpileup" \
        INPUT_MMAP_READS_PILEUP_FILES="../ENCSR238ZZD_ENCFF719MSG_1_ENCFF120MML_2_1_1/ENCFF719MSG_ENCFF120MML_hap1_mmapreads.mpileup ../ENCSR238ZZD_ENCFF719MSG_1_ENCFF120MML_2_1_1/ENCFF719MSG_ENCFF120MML_hap2_mmapreads.mpileup ../ENCSR238ZZD_ENCFF337ZBN_1_ENCFF481IQE_2_1_1/ENCFF337ZBN_ENCFF481IQE_hap1_mmapreads.mpileup ../ENCSR238ZZD_ENCFF337ZBN_1_ENCFF481IQE_2_1_1/ENCFF337ZBN_ENCFF481IQE_hap2_mmapreads.mpileup" \
        PREFIX=ENCSR238ZZD \
        FDR_CUTOFF=0.10 \
        VCF_SAMPLE_ID=sge_Aug_encodev2_2_local
```

#### main output file: ENCSR238ZZD_interestingHets.FDR-0.10.binom.chrs1-22.6-tot_0-min_cnt.tsv



## (3) Aggregate across genomic elements, e.g. genes. ### 

```
make -f ~/bin/AlleleSeq2/PIPELINE_aggregate_over_genomic_regions.mk \
	  PREFIX=ENCSR238ZZD \
	  REGIONS_FILE="../../../../pgenomes_20191127/gencode.v24_all_genes.bed" \
	  COUNTS_FILE=../ENCSR238ZZD_combined/ENCSR238ZZD_filtered_counts.tsv \
	  UNIQ_ALN_FILES='../ENCSR238ZZD_ENCFF719MSG_1_ENCFF120MML_2_1_1/ENCFF719MSG_ENCFF120MML_ASE-params_crdsorted_uniqreads_over_hetSNVs.bam ../ENCSR238ZZD_ENCFF337ZBN_1_ENCFF481IQE_2_1_1/ENCFF337ZBN_ENCFF481IQE_ASE-params_crdsorted_uniqreads_over_hetSNVs.bam' \
	  MMAP_ALN_FILES='../ENCSR238ZZD_ENCFF719MSG_1_ENCFF120MML_2_1_1/ENCFF719MSG_ENCFF120MML_ASE-params_crdsorted_mmapreads_over_hetSNVs.bam ../ENCSR238ZZD_ENCFF337ZBN_1_ENCFF481IQE_2_1_1/ENCFF337ZBN_ENCFF481IQE_ASE-params_crdsorted_mmapreads_over_hetSNVs.bam' \
	  FDR_CUTOFF=0.10  
```


#### all input files should be produced from (1) or (2)

####  REGIONS_FILE -- must be a bed file: e.g., coordinates and gene name

#### main output file with ASE genes: ENCSR238ZZD_interesting_regions.FDR-0.10.betabinom.chrs1-22.6-tot_0-min_cnt.tsv 
