#######################
### PIPELINE PARAMS ###
#######################


### system / executables ##

PL                               := ~/bin/AlleleSeq2
PGENOME_DIR                      := NULL
VCF_SAMPLE_ID                    := NULL
INPUT_UNIQ_READS_PILEUP_FILES    := NULL # list of .mpileup files with uniquely mapped reads to aggregate
INPUT_MMAP_READS_PILEUP_FILES    := NULL # list of .mpileup files with multi-mapping reads
PREFIX                           := NULL

# stats / counts params:

FDR_SIMS                         := 500
FDR_CUTOFF                       := 0.05
Cntthresh                        := 6
AMB_MODE                         := adjust # 'adjust' or allelic ratio diff threshold for filtering



######################
### PIPELINE STEPS ###
######################


$(info PGENOME_DIR: $(PGENOME_DIR))
$(info INPUT_UNIQ_READS_PILEUP_FILES: $(INPUT_UNIQ_READS_PILEUP_FILES))
$(info INPUT_MMAP_READS_PILEUP_FILES: $(INPUT_MMAP_READS_PILEUP_FILES))
$(info PREFIX: $(PREFIX))
$(info $(empty_string))


######################
### PIPELINE START ###
######################

all: $(PREFIX)_raw_counts_ref_allele_ratios.pdf $(PREFIX)_final_counts_ref_allele_ratios.pdf interestingHets.binom.tsv interestingHets.betabinom.tsv 

# this seems to work, but the way it deals with paths, filenames, etc needs to be cleaned up
interestingHets.betabinom.tsv: $(PREFIX)_final_counts.tsv 
	Rscript $(PL)/alleledb_calcOverdispersion.R $< betabinomial
	Rscript $(PL)/alleledb_alleleseqBetabinomial.R $< betabinomial counts.betabinom.tsv interestingHets.betabinom.tsv FDR.betabinomial.txt $(FDR_CUTOFF)

interestingHets.binom.tsv: $(PREFIX)_final_counts.tsv
	python $(PL)/FalsePos.py $< $(FDR_SIMS) $(FDR_CUTOFF) > FDR.binom.txt
	cat $< | python $(PL)/filter_by_pval.py FDR.binom.txt > $@

# allelic ratio distrs
$(PREFIX)_final_counts_ref_allele_ratios.pdf: $(PREFIX)_final_counts.tsv
	cat $< | python $(PL)/plot_AllelicRatio_distribution.py $(PREFIX)_final_counts

# filter out sites in potential cnv regions and in non-autosomal chr; 
# and sites with seemingly misphased/miscalled nearby variants
# filter/adjust sites imbalanced likely due to unaccounted multi-mapping reads 
$(PREFIX)_final_counts.tsv: $(PREFIX)_raw_counts.tsv $(INPUT_MMAP_READS_PILEUP_FILES)
	cat $< | \
	python $(PL)/filter_cnv_sites.py $(PREFIX)_discarded_HetSNVs.tsv $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_rd.tab | \
	python $(PL)/filter_chr.py $(PREFIX)_discarded_HetSNVs.tsv | \
	python $(PL)/filter_phase_warnings.py $(PREFIX)_discarded_HetSNVs.tsv | \
	python $(PL)/filter_sites_w_mmaps.py $(AMB_MODE) $(PREFIX)_discarded_HetSNVs.tsv filter_sites_w_mmaps.log \
		$(INPUT_MMAP_READS_PILEUP_FILES) > $@

# allelic ratio distrs
$(PREFIX)_raw_counts_ref_allele_ratios.pdf: $(PREFIX)_raw_counts.tsv
	cat $< | python $(PL)/plot_AllelicRatio_distribution.py $(PREFIX)_raw_counts

# counts
$(PREFIX)_raw_counts.tsv: $(INPUT_UNIQ_READS_PILEUP_FILES)
	python $(PL)/pileup2counts.py $(Cntthresh) $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h1.bed $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h2.bed \
	$(PREFIX)_discarded_HetSNVs.tsv \
	$(INPUT_UNIQ_READS_PILEUP_FILES) > $@


