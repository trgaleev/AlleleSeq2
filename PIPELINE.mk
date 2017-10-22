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
#JAVA                             := ~/bin/jre1.8.0_77/bin/java
JAVA                             := java
JAVA_MEM                         := 80g
STAR                             := ~/bin/STAR/bin/Linux_x86_64/STAR


### input files / paths ##

READS_R1                          :=
READS_R2                          :=

PGENOME_DIR                      := NULL
VCF_SAMPLE_ID                    := NULL


### mapping params ##

ALIGNMENT_MODE                           := NULL # can be ASE, ASB, custom
RM_DUPLICATE_READS                       := off  # with 'on' duplicate reads will be removed using picard

# needed for all: ASE, ASB, or custom:
GenomeIdx_STAR_diploid                   := $(PGENOME_DIR)/STAR_idx_diploid
STAR_outFilterMismatchNoverReadLmax      := 0.03
STAR_outFilterMatchNminOverLread         := 0.95
STAR_readFIlesCommand                    := zcat # zcat, cat, etc
STAR_limitSjdbInsertNsj                  := 1500000 # star default is 1000000

# needed if ASE
REFGENOME_VERSION             := GRCh37   #GRCh38 or CRCh37
Annotation_diploid            := $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_diploid.gencode.v19.annotation.gtf

ifeq ($(REFGENOME_VERSION), GRCh38)
  Annotation_diploid          := $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_diploid.gencode.v24.annotation.gtf
endif

STAR_sjdbOverhang                        := 100  #STAR default will work as well as the ideal (readlength -1) value according to the manual

# if custom alignment mode
STAR_parameters_file                     := $(PL)/STAR_custom_parameters_sample_file


# stats / counts params:

FDR_SIMS                         := 500
FDR_CUTOFF                       := 0.05
Cntthresh_total                  := 6
Cntthresh_lower                  := 0
AMB_MODE                         := adjust # 'adjust' or allelic ratio diff threshold for filtering



######################
### PIPELINE STEPS ###
######################

empty_string:=
ifeq ($(READS_R2),$(empty_string))
  tmp1 = $(READS_R1:.gz=)
  tmp2 = $(tmp1:.fastq=)
  PREFIX = $(tmp2:.fq=)
else
  tmp11 = $(READS_R1:.gz=)
  tmp12 = $(tmp11:.fastq=)
  tmp21 = $(READS_R2:.gz=)
  tmp22 = $(tmp21:.fastq=)
  PREFIX = $(tmp12:.fq=)_$(tmp22:.fq=)
endif

ifeq ($(RM_DUPLICATE_READS),on)
  DEDUP_SUFFIX = rmdup.
endif


FINAL_ALIGNMENT_FILENAME = $(PREFIX)_$(ALIGNMENT_MODE)-params.Aligned.sortedByCoord.out.$(DEDUP_SUFFIX)bam
HetSNV_UNIQALNS_FILENAME = $(PREFIX)_$(ALIGNMENT_MODE)-params_crdsorted_uniqreads_over_hetSNVs.bam
HetSNV_MMAPALNS_FILENAME = $(PREFIX)_$(ALIGNMENT_MODE)-params_crdsorted_mmapreads_over_hetSNVs.bam

$(info PGENOME_DIR: $(PGENOME_DIR))
$(info READS_R1: $(READS_R1))
$(info READS_R2: $(READS_R2))
$(info PREFIX: $(PREFIX))
$(info ALIGNMENT_MODE: $(ALIGNMENT_MODE))
$(info RM_DUPLICATE_READS: $(RM_DUPLICATE_READS))
$(info FINAL_ALIGNMENT_FILENAME: $(FINAL_ALIGNMENT_FILENAME))
$(info HetSNV_UNIQALNS_FILENAME: $(HetSNV_UNIQALNS_FILENAME))
$(info HetSNV_MMAPALNS_FILENAME: $(HetSNV_MMAPALNS_FILENAME))
$(info $(empty_string))



######################
### PIPELINE START ###
######################

all: $(PREFIX)_raw_counts_ref_allele_ratios.pdf $(PREFIX)_final_counts_ref_allele_ratios.pdf interestingHets.binom.tsv interestingHets.betabinom.tsv

# this seems to work, but the way it deals with paths, filenames, etc needs to be cleaned up
interestingHets.betabinom.tsv: $(PREFIX)_final_counts_min.$(Cntthresh_total)-total.$(Cntthresh_lower)-lower.tsv 
	Rscript $(PL)/alleledb_calcOverdispersion.R $< betabinomial
	Rscript $(PL)/alleledb_alleleseqBetabinomial.R $< betabinomial counts.betabinom.tsv interestingHets.betabinom.tsv FDR.betabinomial.txt $(FDR_CUTOFF)

interestingHets.binom.tsv: $(PREFIX)_final_counts_min.$(Cntthresh_total)-total.$(Cntthresh_lower)-lower.tsv
	python $(PL)/FalsePos.py $< $(FDR_SIMS) $(FDR_CUTOFF) > FDR.binom.txt
	cat $< | python $(PL)/filter_by_pval.py FDR.binom.txt > $@


# allelic ratio distrs
$(PREFIX)_final_counts_min.$(Cntthresh_total)-total.$(Cntthresh_lower)-lower_ref_allele_ratios.pdf: $(PREFIX)_final_counts_min.$(Cntthresh_total)-total.$(Cntthresh_lower)-lower.tsv
	cat $< | python $(PL)/plot_AllelicRatio_distribution.py $(PREFIX)_final_counts_min.$(Cntthresh_total)-total.$(Cntthresh_lower)-lower

$(PREFIX)_final_counts_min.$(Cntthresh_total)-total.$(Cntthresh_lower)-lower.tsv: $(PREFIX)_final_counts.tsv
	cat $< | python $(PL)/filter_by_counts.py $(Cntthresh_total) $(Cntthresh_lower) > $@


# allelic ratio distrs
$(PREFIX)_final_counts_ref_allele_ratios.pdf: $(PREFIX)_final_counts.tsv
	cat $< | python $(PL)/plot_AllelicRatio_distribution.py $(PREFIX)_final_counts 
# filter out sites in potential cnv regions and in non-autosomal chr; 
# and sites with seemingly misphased/miscalled nearby variants
# filter/adjust sites imbalanced likely due to unaccounted multi-mapping reads 
$(PREFIX)_final_counts.tsv: $(PREFIX)_raw_counts.tsv $(PREFIX)_h1_mmapreads.mpileup $(PREFIX)_h2_mmapreads.mpileup
	cat $< | \
	python $(PL)/filter_cnv_sites.py $(PREFIX)_discarded_HetSNVs.tsv $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_rd.tab | \
	python $(PL)/filter_chr.py $(PREFIX)_discarded_HetSNVs.tsv | \
	python $(PL)/filter_phase_warnings.py $(PREFIX)_discarded_HetSNVs.tsv | \
	python $(PL)/filter_sites_w_mmaps.py $(AMB_MODE) $(PREFIX)_discarded_HetSNVs.tsv filter_sites_w_mmaps.log \
		$(PREFIX)_h1_mmapreads.mpileup $(PREFIX)_h2_mmapreads.mpileup > $@

# allelic ratio distrs
$(PREFIX)_raw_counts_ref_allele_ratios.pdf: $(PREFIX)_raw_counts.tsv
	cat $< | python $(PL)/plot_AllelicRatio_distribution.py $(PREFIX)_raw_counts

# counts
$(PREFIX)_raw_counts.tsv: $(PREFIX)_h1_uniqreads.mpileup $(PREFIX)_h2_uniqreads.mpileup
	python $(PL)/pileup2counts.py 1 $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h1.bed $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h2.bed \
	$(PREFIX)_discarded_HetSNVs.tsv \
	$(PREFIX)_h1_uniqreads.mpileup $(PREFIX)_h2_uniqreads.mpileup > $@


# pileups
$(PREFIX)_h1_mmapreads.mpileup: $(HetSNV_MMAPALNS_FILENAME)
	$(SAMTOOLS) mpileup -BQ0 --max-depth 999999 --ff UNMAP -f $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_h1.fa $< \
	--positions $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h1.bed --output $@

$(PREFIX)_h2_mmapreads.mpileup: $(HetSNV_MMAPALNS_FILENAME)
	$(SAMTOOLS) mpileup -BQ0 --max-depth 999999 --ff UNMAP -f $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_h2.fa $< \
	--positions $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h2.bed --output $@

$(PREFIX)_h1_uniqreads.mpileup: $(HetSNV_UNIQALNS_FILENAME)
	$(SAMTOOLS) mpileup -BQ0 --max-depth 999999 --ff UNMAP -f $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_h1.fa $< \
	--positions $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h1.bed --output $@

$(PREFIX)_h2_uniqreads.mpileup: $(HetSNV_UNIQALNS_FILENAME)
	$(SAMTOOLS) mpileup -BQ0 --max-depth 999999 --ff UNMAP -f $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_h2.fa $< \
	--positions $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h2.bed --output $@

# by default --ff was also filtering-out some other - secondary? reads:  [UNMAP,SECONDARY,QCFAIL,DUP], leaving only UNMAP for now




# non-uniq alns over hetSNVs:
$(PREFIX)_$(ALIGNMENT_MODE)-params_crdsorted_mmapreads_over_hetSNVs.bam: $(FINAL_ALIGNMENT_FILENAME)
	cat $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h1.bed $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h2.bed | \
	$(SAMTOOLS) view -h -L - $< | awk '$$5!="255" {print $$0}' | \
	$(SAMTOOLS) view -b - > $@
	$(SAMTOOLS) index $@
	$(SAMTOOLS) flagstat $@ > $@.stat

# uniq alns over hetSNVs:
$(PREFIX)_$(ALIGNMENT_MODE)-params_crdsorted_uniqreads_over_hetSNVs.bam: $(FINAL_ALIGNMENT_FILENAME)
	cat $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h1.bed $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h2.bed | \
	$(SAMTOOLS) view -h -q 255 -L - $< | \
	$(SAMTOOLS) view -b - > $@
	$(SAMTOOLS) index $@
	$(SAMTOOLS) flagstat $@ > $@.stat


# if removing duplicate reads
$(PREFIX)_$(ALIGNMENT_MODE)-params.Aligned.sortedByCoord.out.rmdup.bam: $(PREFIX)_$(ALIGNMENT_MODE)-params.Aligned.sortedByCoord.out.bam
	$(JAVA) -Xmx$(JAVA_MEM) -jar $(PICARD) MarkDuplicates \
	INPUT=$< OUTPUT=$@ METRICS_FILE=$(@:.bam=.metrics) \
	REMOVE_DUPLICATES=true \
	DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES
	$(SAMTOOLS) index $@
	$(SAMTOOLS) flagstat $@ > $@.stat


## specific additional params will be read from $(STAR_parameters_file)
$(PREFIX)_custom-params.Aligned.sortedByCoord.out.bam: $(READS_R1)
	$(STAR) \
	--runThreadN $(NTHR) \
	--genomeDir $(GenomeIdx_STAR_diploid) \
	--readFilesIn $< $(READS_R2) \
	--readFilesCommand $(STAR_readFIlesCommand) \
	--outFileNamePrefix $(@:Aligned.sortedByCoord.out.bam=) \
	--outSAMattributes All \
	--outFilterMultimapNmax 999999 \
	--scoreGenomicLengthLog2scale 0.0 \
	--sjdbOverhang $(STAR_sjdbOverhang) \
	--sjdbGTFfile $(Annotation_diploid) \
	--parametersFiles $(STAR_parameters_file) \
	--limitSjdbInsertNsj $(STAR_limitSjdbInsertNsj) \
	--outSAMtype BAM SortedByCoordinate
	$(SAMTOOLS) flagstat $@ > $@.stat
	$(SAMTOOLS) index $@

## optimal? params for RNA-seq; will use as default for ASE
$(PREFIX)_ASE-params.Aligned.sortedByCoord.out.bam: $(READS_R1)
	$(STAR) \
	--runThreadN $(NTHR) \
	--genomeDir $(GenomeIdx_STAR_diploid) \
	--twopassMode Basic \
	--readFilesIn $< $(READS_R2) \
	--readFilesCommand $(STAR_readFIlesCommand) \
	--outFileNamePrefix $(@:Aligned.sortedByCoord.out.bam=) \
	--outSAMattributes All \
	--outFilterMismatchNoverReadLmax $(STAR_outFilterMismatchNoverReadLmax) \
	--outFilterMatchNminOverLread $(STAR_outFilterMatchNminOverLread) \
	--outFilterMultimapNmax 999999 \
	--scoreGenomicLengthLog2scale 0.0 \
	--sjdbOverhang $(STAR_sjdbOverhang) \
	--sjdbGTFfile $(Annotation_diploid) \
	--limitSjdbInsertNsj $(STAR_limitSjdbInsertNsj) \
	--outSAMtype BAM SortedByCoordinate
	$(SAMTOOLS) flagstat $@ > $@.stat
	$(SAMTOOLS) index $@	

## opts similar to AlleleSeq v1.2a bowtie1 -v 2 -m 1 mode; but with small gaps allowed, no splicing; will use as default for ASB
$(PREFIX)_ASB-params.Aligned.sortedByCoord.out.bam: $(READS_R1)
	$(STAR) \
	--runThreadN $(NTHR) \
	--genomeDir $(GenomeIdx_STAR_diploid) \
	--readFilesIn $< $(READS_R2) \
	--readFilesCommand $(STAR_readFIlesCommand) \
	--outFileNamePrefix $(@:Aligned.sortedByCoord.out.bam=) \
	--outSAMattributes All \
	--outFilterMismatchNoverReadLmax $(STAR_outFilterMismatchNoverReadLmax) \
	--outFilterMatchNminOverLread $(STAR_outFilterMatchNminOverLread) \
	--outFilterMultimapNmax 999999 \
	--scoreGap -100 \
	--scoreGenomicLengthLog2scale 0.0 \
	--sjdbScore 0 \
	--limitSjdbInsertNsj $(STAR_limitSjdbInsertNsj) \
	--outSAMtype BAM SortedByCoordinate
	$(SAMTOOLS) flagstat $@ > $@.stat
	$(SAMTOOLS) index $@	



# * --outFilterMultimapScoreRange default =1 works fine for making sure read mapping to the right allele with one hetSNV gets a higher score
# 	and treated as uniquely mapped: the score diff between a perfectly mapped read and one with one mismatch is 2 (with --scoreGenomicLengthLog2scale 0.0):
# 	one less mapped base and -1 for mismatch?
#



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
