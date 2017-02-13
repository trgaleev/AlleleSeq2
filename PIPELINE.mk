#######################
### PIPELINE PARAMS ###
#######################

PL                               := ~/bin/AlleleSeq2
BOWTIE1                          := ~/bin/bowtie-1.1.1/bowtie
NTHR                             := 1 # multithread, works for mapping, sorting
BOWTIE1                          := /home/fas/gerstein/tg397/distr/bowtie-1.1.1/bowtie
SAMTOOLS                         := ~/bin/samtools-1.3.1/samtools
PICARD                           := ~/distr/picard-tools-2.1.1/picard.jar
JAVA                             := ~/distr/jre1.8.0_77/bin/java
JAVA_MEM                         := 80g
STAR                             := ~/bin/STAR/bin/Linux_x86_64/STAR
READS_R1                         := NULL
READS_R2                         := NULL
READS                            := NULL
PGENOME_DIR                      := NULL
VCF_SAMPLE_ID                    := NULL
STAR_sjdbOverhang                    := 100
STAR_outFilterMismatchNoverReadLmax  := 0.03
STAR_outFilterMatchNminOverLread     := 0.95
STAR_parameters_file                 := $(PL)/STAR_custom_parameters_sample_file
STAR_readFIlesCommand            := zcat # zcat, cat, etc
ALIGNMENT_MODE                   := NULL # can be bowtie1, STAR_mode-0, STAR_mode-1, STAR_mode-2, STAR_mode-custom
AMBIGOUS_MAPPING_BIAS            := filter_out # ignore, filter_out, adjust_counts
FDR_SIMS                         := 5
FDR_CUTOFF                       := 0.05
Cntthresh                        := 6


######################
### PIPELINE STEPS ###
######################

ifeq ($(READS_R2),NULL)
  tmp = $(READS:.fq.gz=)
  PREFIX = $(tmp:.fastq.gz=)
else
  PREFIX = _temp
  SUFFIX_PE = _R1R2
endif

ifeq ($(AMBIGOUS_MAPPING_BIAS),ignore)
  SUFFIX_AMB = _uniqonly
else ifeq ($(SUFFIX_AMB), filter_out)
  SUFFIX_AMB = _temp1
else ifeq ($(SUFFIX_AMB), adjust_counts)
  SUFFIX_AMB = _temp2
endif

### something strange happens with bowtie mapping to the large genome, so turning this off for now! 
#ifeq ($(ALIGNMENT_MODE),bowtie1)
#  FINAL_BAM_FILE = $(PREFIX)_$(ALIGNMENT_MODE)$(SUFFIX_PE)$(SUFFIX_AMB)_coordsorted
#  GenomeIdx_bowtie_diploid = $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_diploid
#else
#  FINAL_BAM_FILE = $(PREFIX)_$(ALIGNMENT_MODE)$(SUFFIX_PE)$(SUFFIX_AMB).Aligned.sortedByCoord.out
#  GenomeIdx_STAR_diploid = $(PGENOME_DIR)/STAR_idx_diploid
#endif

FINAL_BAM_FILE = $(PREFIX)_$(ALIGNMENT_MODE)$(SUFFIX_PE)$(SUFFIX_AMB).Aligned.sortedByCoord.out
GenomeIdx_STAR_diploid := $(PGENOME_DIR)/STAR_idx_diploid
Annotation_diploid := $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_diploid.gencode.v19.annotation.gtf
FINAL_BAM_FILE := $(FINAL_BAM_FILE).bam



$(warning READS: $(READS))
$(warning FINAL_BAM_FILE: $(FINAL_BAM_FILE))
$(warning ALIGNMENT_MODE: $(ALIGNMENT_MODE))
$(warning PREFIX: $(PREFIX))
$(warning GenomeIdx_bowtie_diploid: $(GenomeIdx_bowtie_diploid))
$(warning SUFFIX_PE: $(SUFFIX_PE))
$(warning AMBIGOUS_MAPPING_BIAS: $(AMBIGOUS_MAPPING_BIAS))
$(warning SUFFIX_AMB: $(SUFFIX_AMB))

######################
### PIPELINE START ###
######################

all: $(PREFIX)_counts_noCNV_noXYM.tsv $(PREFIX)_counts_ref_allele_ratios.pdf  $(PREFIX)_counts_noCNV_noXYM_ref_allele_ratios.pdf interestingHets.betabinom.tsv

# this seems to work, but the way it deals with paths, filenames, etc needs to be cleaned up
interestingHets.betabinom.tsv: $(PREFIX)_counts_noCNV_noXYM.tsv
	cat $< | Rscript $(PL)/alleledb_calcOverdispersion.R $(CURDIR) $(PREFIX)_counts_noCNV_noXYM.tsv
	cat $< | Rscript $(PL)/alleledb_alleleseqBetabinomial.R $(CURDIR) $(PREFIX)_counts_noCNV_noXYM.tsv $(FDR_CUTOFF)

# allelic ratio distrs
$(PREFIX)_counts_noCNV_noXYM_ref_allele_ratios.pdf: $(PREFIX)_counts_noCNV_noXYM.tsv
	cat $< | Rscript $(PL)/counts_allelic_ratio_distribution_plot.R $(PREFIX)_filtered_counts

# keep autosomal chr only
$(PREFIX)_counts_noCNV_noXYM.tsv: $(PREFIX)_counts_noCNV.tsv
	cat $< | python $(PL)/filter_chr.py > $@

# rm CNV-sites
$(PREFIX)_counts_noCNV.tsv: $(PREFIX)_counts.tsv
	python $(PL)/filter_cnv_sites.py $(PGENOME_DIR)/$(VCF_SAMPLE_ID).alleleSeqInput.cnv $< > $@

#allelic ratio distrs
$(PREFIX)_counts_ref_allele_ratios.pdf: $(PREFIX)_counts.tsv
	cat $< | Rscript $(PL)/counts_allelic_ratio_distribution_plot.R $(PREFIX)_counts

# counts
$(PREFIX)_counts.tsv: $(PREFIX)_h1.mpileup $(PREFIX)_h2.mpileup
	python $(PL)/mpileup2counts.py $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h1.bed $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h2.bed $(PREFIX)_h1.mpileup $(PREFIX)_h2.mpileup $(Cntthresh) > $@

# pileups
$(PREFIX)_h1.mpileup: $(FINAL_BAM_FILE)
	$(SAMTOOLS) mpileup -BQ0 --max-depth 999999 -f $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_h1.fa $< --positions $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h1.bed --output $@

$(PREFIX)_h2.mpileup: $(FINAL_BAM_FILE)
	$(SAMTOOLS) mpileup -BQ0 --max-depth 999999 -f $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_h2.fa $< --positions $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h2.bed --output $@


$(PREFIX)_STAR_mode-custom_uniqonly.Aligned.sortedByCoord.out.bam: $(READS)
	$(STAR) \
	--runThreadN $(NTHR) \
	--genomeDir $(GenomeIdx_STAR_diploid) \
	--readFilesIn $< \
	--readFilesCommand $(STAR_readFIlesCommand) \
	--outFileNamePrefix $(@:Aligned.sortedByCoord.out.bam=) \
	--outSAMattributes All \
	--outFilterMultimapNmax 1 \
	--sjdbOverhang $(STAR_sjdbOverhang) \
	--sjdbGTFfile $(Annotation_diploid) \
	--parametersFiles $(STAR_parameters_file) \
	--outSAMtype BAM SortedByCoordinate
	$(SAMTOOLS) flagstat $@ > $@.stat

$(PREFIX)_STAR_mode-2_uniqonly.Aligned.sortedByCoord.out.bam: $(READS)
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
	--outFilterMultimapNmax 1 \
	--sjdbOverhang $(STAR_sjdbOverhang) \
	--sjdbGTFfile $(Annotation_diploid) \
	--outSAMtype BAM SortedByCoordinate
	$(SAMTOOLS) flagstat $@ > $@.stat

$(PREFIX)_STAR_mode-1_uniqonly.Aligned.sortedByCoord.out.bam: $(READS)
	$(STAR) \
	--runThreadN $(NTHR) \
	--genomeDir $(GenomeIdx_STAR_diploid) \
	--readFilesIn $< \
	--readFilesCommand $(STAR_readFIlesCommand) \
	--outFileNamePrefix $(@:Aligned.sortedByCoord.out.bam=) \
	--outSAMattributes All \
	--outFilterMismatchNoverReadLmax $(STAR_outFilterMismatchNoverReadLmax) \
	--outFilterMatchNminOverLread $(STAR_outFilterMatchNminOverLread) \
	--outFilterMultimapNmax 1 \
	--scoreGap -100 \
	--scoreGenomicLengthLog2scale 0.0 \
	--sjdbScore 0 \
	--outSAMtype BAM SortedByCoordinate
	$(SAMTOOLS) flagstat $@ > $@.stat

$(PREFIX)_STAR_mode-0_uniqonly.Aligned.sortedByCoord.out.bam: $(READS)
	$(STAR) \
	--runThreadN $(NTHR) \
	--genomeDir $(GenomeIdx_STAR_diploid) \
	--readFilesIn $< \
	--readFilesCommand $(STAR_readFIlesCommand) \
	--outFileNamePrefix $(@:Aligned.sortedByCoord.out.bam=) \
	--outSAMattributes All \
	--outFilterMultimapNmax 1 \
	--outFilterMismatchNmax 2 \
	--alignEndsType EndToEnd \
	--scoreDelOpen -100 \
	--scoreInsOpen -100 \
	--scoreGap -100 \
	--scoreGenomicLengthLog2scale 0.0 \
	--sjdbScore 0 \
	--outSAMtype BAM SortedByCoordinate
	$(SAMTOOLS) flagstat $@ > $@.stat

$(PREFIX)_bowtie1_uniqonly_coordsorted.bam: $(READS)
	zcat $(READS) | \
	$(BOWTIE1) -S -p $(NTHR) --best --strata -v 2 -m 1 $(GenomeIdx_bowtie_diploid) --large-index - | \
	$(SAMTOOLS) view -F 0x04 -b - | \
	$(SAMTOOLS) sort - -o $@
	$(SAMTOOLS) flagstat $@ > $@.stat



