#######################
### PIPELINE PARAMS ###
#######################


### system / executables ##

PL                               := ~/bin/AlleleSeq2
#BOWTIE1                          := ~/bin/bowtie-1.1.1/bowtie
NTHR                             := 1 # multithread, works for mapping, sorting
#BOWTIE1                          := /home/fas/gerstein/tg397/distr/bowtie-1.1.1/bowtie
SAMTOOLS                         := ~/bin/samtools-1.3.1/samtools
PICARD                           := ~/bin/picard-tools-2.1.1/picard.jar
JAVA                             := ~/bin/jre1.8.0_77/bin/java
JAVA_MEM                         := 80g
STAR                             := ~/bin/STAR/bin/Linux_x86_64/STAR


### input files / paths ##

READS_1                          := NULL 
READS_2                          := NULL
READS                            := NULL # specify READS if SE, or READS_1 and READS_2 if PE (the latter doesn't work yet)

PGENOME_DIR                      := NULL
VCF_SAMPLE_ID                    := NULL


### mapping params ##

ALIGNMENT_MODE                           := NULL # can be ASE, ASB, custom
RM_DUPs                                  := off  # 'on' doesn't work yet

# needed for all: ASE, ASB, or custom:
GenomeIdx_STAR_diploid                   := $(PGENOME_DIR)/STAR_idx_diploid
STAR_outFilterMismatchNoverReadLmax      := 0.03
STAR_outFilterMatchNminOverLread         := 0.95
STAR_readFIlesCommand                    := zcat # zcat, cat, etc

# needed if ASE
Annotation_diploid                       := $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_diploid.gencode.v19.annotation.gtf
STAR_sjdbOverhang                        := 100  #STAR default will work as well as the ideal (readlength -1) value according to the manual

# if custom alignment mode
STAR_parameters_file                     := $(PL)/STAR_custom_parameters_sample_file


# stats / counts params:

FDR_SIMS                         := 500
FDR_CUTOFF                       := 0.05
Cntthresh                        := 6
AMB_MODE                         := adjust # 'adjust' or allelic ratio diff threshold for filtering

######################
### PIPELINE STEPS ###
######################

ifeq ($(READS_2),NULL)
  tmp = $(READS:.fq.gz=)
  PREFIX = $(tmp:.fastq.gz=)
else
  tmp = $(READS_2:.fq.gz=)
  tmp1 = $(tmp:.fastq.gz=)
  PREFIX = $(tmp1:_2=)
endif


ALIGMENT_FILENAME         = $(PREFIX)_$(ALIGNMENT_MODE)-params.Aligned.sortedByCoord.out.bam
HetSNV_UNIQALNS_FILENAME  = $(PREFIX)_$(ALIGNMENT_MODE)-params_crdsorted_uniqreads_over_hetSNVs.bam
HetSNV_MMAPALNS_FILENAME  = $(PREFIX)_$(ALIGNMENT_MODE)-params_crdsorted_mmapreads_over_hetSNVs.bam

$(info READS: $(READS))
$(info ALIGMENT_FILENAME: $(ALIGMENT_FILENAME))
$(info ALIGNMENT_MODE: $(ALIGNMENT_MODE))
$(info PREFIX: $(PREFIX))

######################
### PIPELINE START ###
######################

all: $(PREFIX)_raw_counts_ref_allele_ratios.pdf $(PREFIX)_final_counts_ref_allele_ratios.pdf interestingHets.betabinom.tsv

# this seems to work, but the way it deals with paths, filenames, etc needs to be cleaned up
interestingHets.betabinom.tsv: $(PREFIX)_final_counts.tsv 
	cat $< | Rscript $(PL)/alleledb_calcOverdispersion.R $(CURDIR) $(PREFIX)_final_counts.tsv
	cat $< | Rscript $(PL)/alleledb_alleleseqBetabinomial.R $(CURDIR) $(PREFIX)_final_counts.tsv $(FDR_CUTOFF)

# allelic ratio distrs
$(PREFIX)_final_counts_ref_allele_ratios.pdf: $(PREFIX)_final_counts.tsv
	cat $< | Rscript $(PL)/counts_allelic_ratio_distribution_plot.R $(PREFIX)_final_counts

# filter out sites in potential cnv regions and in non-autosomal chr; 
# filter/adjust sites imbalanced likely due to unaccounted multi-mapping reads 
# and sites with seemingly misphased/miscalled nearby variants
$(PREFIX)_final_counts.tsv: $(PREFIX)_raw_counts.tsv $(PREFIX)_h1_mmapreads.mpileup $(PREFIX)_h2_mmapreads.mpileup
	cat $< | \
	python $(PL)/filter_cnv_sites.py $(PREFIX)_discarded_HetSNVs.tsv $(PGENOME_DIR)/$(VCF_SAMPLE_ID).alleleSeqInput.cnv | \
	python $(PL)/filter_chr.py $(PREFIX)_discarded_HetSNVs.tsv | \
	python $(PL)/filter_phase_warnings.py $(PREFIX)_discarded_HetSNVs.tsv | \
	python $(PL)/filter_sites_w_mmaps.py $(PREFIX)_h1_mmapreads.mpileup $(PREFIX)_h2_mmapreads.mpileup $(AMB_MODE) $(PREFIX)_discarded_HetSNVs.tsv filter_mm.log > $@

# allelic ratio distrs
$(PREFIX)_raw_counts_ref_allele_ratios.pdf: $(PREFIX)_raw_counts.tsv
	cat $< | Rscript $(PL)/counts_allelic_ratio_distribution_plot.R $(PREFIX)_raw_counts

# counts
$(PREFIX)_raw_counts.tsv: $(PREFIX)_h1_uniqreads.mpileup $(PREFIX)_h2_uniqreads.mpileup
	python $(PL)/pileup2counts.py $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h1.bed $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h2.bed $(PREFIX)_h1_uniqreads.mpileup $(PREFIX)_h2_uniqreads.mpileup $(Cntthresh) > $@


# pileups
$(PREFIX)_h1_mmapreads.mpileup: $(HetSNV_MMAPALNS_FILENAME)
	$(SAMTOOLS) mpileup -BQ0 --max-depth 999999 --ff UNMAP -f $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_h1.fa $< --positions $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h1.bed --output $@

$(PREFIX)_h2_mmapreads.mpileup: $(HetSNV_MMAPALNS_FILENAME)
	$(SAMTOOLS) mpileup -BQ0 --max-depth 999999 --ff UNMAP -f $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_h2.fa $< --positions $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h2.bed --output $@

$(PREFIX)_h1_uniqreads.mpileup: $(HetSNV_UNIQALNS_FILENAME)
	$(SAMTOOLS) mpileup -BQ0 --max-depth 999999 --ff UNMAP -f $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_h1.fa $< --positions $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h1.bed --output $@

$(PREFIX)_h2_uniqreads.mpileup: $(HetSNV_UNIQALNS_FILENAME)
	$(SAMTOOLS) mpileup -BQ0 --max-depth 999999 --ff UNMAP -f $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_h2.fa $< --positions $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h2.bed --output $@

# by default --ff was also filtering-out some other - secondary? reads:  [UNMAP,SECONDARY,QCFAIL,DUP], leaving only UNMAP for now


# non-uniq alns over hetSNVs:
$(HetSNV_MMAPALNS_FILENAME): $(ALIGMENT_FILENAME)
	cat $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h1.bed $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h2.bed | \
	$(SAMTOOLS) view -h -L - $< | grep -wv 'NH:i:1' | \
	$(SAMTOOLS) view -b - > $@
	$(SAMTOOLS) index $@

# uniq alns over hetSNVs:
$(HetSNV_UNIQALNS_FILENAME): $(ALIGMENT_FILENAME)
	cat $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h1.bed $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h2.bed | \
	$(SAMTOOLS) view -h -L - $< | grep -w '^@..\|NH:i:1' | \
	$(SAMTOOLS) view -b - > $@
	$(SAMTOOLS) index $@

## specific additional params will be read from $(STAR_parameters_file)
$(PREFIX)_custom-params.Aligned.sortedByCoord.out.bam: $(READS)
	$(STAR) \
	--runThreadN $(NTHR) \
	--genomeDir $(GenomeIdx_STAR_diploid) \
	--readFilesIn $< \
	--readFilesCommand $(STAR_readFIlesCommand) \
	--outFileNamePrefix $(@:Aligned.sortedByCoord.out.bam=) \
	--outSAMattributes All \
	--outFilterMultimapNmax 999999 \
	--sjdbOverhang $(STAR_sjdbOverhang) \
	--sjdbGTFfile $(Annotation_diploid) \
	--parametersFiles $(STAR_parameters_file) \
	--outSAMtype BAM SortedByCoordinate
	$(SAMTOOLS) flagstat $@ > $@.stat
	$(SAMTOOLS) index $@

## optimal? params for RNA-seq; will use as default for ASE
$(PREFIX)_ASE-params.Aligned.sortedByCoord.out.bam: $(READS)
	$(STAR) \
	--runThreadN $(NTHR) \
	--genomeDir $(GenomeIdx_STAR_diploid) \
	--twopassMode Basic \
	--readFilesIn $< \
	--readFilesCommand $(STAR_readFIlesCommand) \
	--outFileNamePrefix $(@:Aligned.sortedByCoord.out.bam=) \
	--outSAMattributes All \
	--outFilterMismatchNoverReadLmax $(STAR_outFilterMismatchNoverReadLmax) \
	--outFilterMatchNminOverLread $(STAR_outFilterMatchNminOverLread) \
	--outFilterMultimapNmax 999999 \
	--sjdbOverhang $(STAR_sjdbOverhang) \
	--sjdbGTFfile $(Annotation_diploid) \
	--outSAMtype BAM SortedByCoordinate
	$(SAMTOOLS) flagstat $@ > $@.stat
	$(SAMTOOLS) index $@	

## opts similar to AlleleSeq v1.2a bowtie1 -v 2 -m 1 mode; but with small gaps allowed, no splicing; will use as default for ASB
$(PREFIX)_ASB-params.Aligned.sortedByCoord.out.bam: $(READS)
	$(STAR) \
	--runThreadN $(NTHR) \
	--genomeDir $(GenomeIdx_STAR_diploid) \
	--readFilesIn $< \
	--readFilesCommand $(STAR_readFIlesCommand) \
	--outFileNamePrefix $(@:Aligned.sortedByCoord.out.bam=) \
	--outSAMattributes All \
	--outFilterMismatchNoverReadLmax $(STAR_outFilterMismatchNoverReadLmax) \
	--outFilterMatchNminOverLread $(STAR_outFilterMatchNminOverLread) \
	--outFilterMultimapNmax 999999 \
	--scoreGap -100 \
	--scoreGenomicLengthLog2scale 0.0 \
	--sjdbScore 0 \
	--outSAMtype BAM SortedByCoordinate
	$(SAMTOOLS) flagstat $@ > $@.stat
	$(SAMTOOLS) index $@	


## opts similar to AlleleSeq v1.2a bowtie1 -v 2 -m 1 mode
#$(PREFIX)_STAR_mode-0_uniqonly.Aligned.sortedByCoord.out.bam: $(READS)
#	$(STAR) \
#	--runThreadN $(NTHR) \
#	--genomeDir $(GenomeIdx_STAR_diploid) \
#	--readFilesIn $< \
#	--readFilesCommand $(STAR_readFIlesCommand) \
#	--outFileNamePrefix $(@:Aligned.sortedByCoord.out.bam=) \
#	--outSAMattributes All \
#	--outFilterMultimapNmax 1 \
#	--outFilterMismatchNmax 2 \
#	--alignEndsType EndToEnd \
#	--scoreDelOpen -100 \
#	--scoreInsOpen -100 \
#	--scoreGap -100 \
#	--scoreGenomicLengthLog2scale 0.0 \
#	--sjdbScore 0 \
#	--outSAMtype BAM SortedByCoordinate
#	$(SAMTOOLS) flagstat $@ > $@.stat
#	$(SAMTOOLS) index $@	

## this doesn't seem to work for some reason - observing strange read alignments with the large index
#$(PREFIX)_bowtie1_uniqonly_coordsorted.bam: $(READS)
#	zcat $(READS) | \
#	$(BOWTIE1) -S -p $(NTHR) --best --strata -v 2 -m 1 $(GenomeIdx_bowtie_diploid) --large-index - | \
#	$(SAMTOOLS) view -F 0x04 -b - | \
#	$(SAMTOOLS) sort - -o $@
#	$(SAMTOOLS) flagstat $@ > $@.stat
