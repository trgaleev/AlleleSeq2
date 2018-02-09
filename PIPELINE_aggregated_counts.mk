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
Cntthresh_tot                  := 6
Cntthresh_min                  := 0
AMB_MODE                         := adjust # 'adjust' or allelic ratio diff threshold for filtering
KEEP_CHR                         := # empty or 'X'



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

all: $(PREFIX)_ref_allele_ratios.raw_counts.pdf $(PREFIX)_ref_allele_ratios.filtered_counts.chrs1-22$(KEEP_CHR).pdf $(PREFIX)_interestingHets.FDR-$(FDR_CUTOFF).binom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv $(PREFIX)_interestingHets.FDR-$(FDR_CUTOFF).betabinom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv 

#all: $(PREFIX)_interestingHets.FDR-$(FDR_CUTOFF).betabinom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv 

# this seems to work, but the way it deals with paths, filenames, etc needs to be cleaned up
$(PREFIX)_interestingHets.FDR-$(FDR_CUTOFF).betabinom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv: $(PREFIX)_filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv
	Rscript $(PL)/alleledb_calcOverdispersion.R $< $(PREFIX)_FDR-$(FDR_CUTOFF).betabinomial.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min 
	Rscript $(PL)/alleledb_alleleseqBetabinomial.R $< $(PREFIX)_FDR-$(FDR_CUTOFF).betabinomial.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min $(PREFIX)_counts.FDR-$(FDR_CUTOFF).betabinom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv $@ $(PREFIX)_FDR-$(FDR_CUTOFF).betabinom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.txt $(FDR_CUTOFF)

$(PREFIX)_interestingHets.FDR-$(FDR_CUTOFF).binom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv: $(PREFIX)_filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv
	python $(PL)/FalsePos.py $< $(FDR_SIMS) $(FDR_CUTOFF) > $(PREFIX)_FDR-$(FDR_CUTOFF).binom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.txt
	cat $< | python $(PL)/filter_by_pval.py $(PREFIX)_FDR-$(FDR_CUTOFF).binom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.txt > $@


# allelic ratio distrs
$(PREFIX)_ref_allele_ratios.filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min.pdf: $(PREFIX)_filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv
	cat $< | python $(PL)/plot_AllelicRatio_distribution.py $(PREFIX)_filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt

# filter based on total counts and min per allele count
$(PREFIX)_filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv: $(PREFIX)_filtered_counts.chrs1-22$(KEEP_CHR).tsv
	cat $< | python $(PL)/filter_by_counts.py $(Cntthresh_tot) $(Cntthresh_min) > $@
                 
# allelic ratio distrs
$(PREFIX)_ref_allele_ratios.filtered_counts.chrs1-22$(KEEP_CHR).pdf: $(PREFIX)_filtered_counts.chrs1-22$(KEEP_CHR).tsv
	cat $< | python $(PL)/plot_AllelicRatio_distribution.py $(PREFIX)_filtered_counts.chrs1-22$(KEEP_CHR)


# filter out sites in potential cnv regions and 
# in non-autosomal chr; 
# and sites with seemingly misphased/miscalled nearby variants
# filter/adjust sites imbalanced likely due to unaccounted multi-mapping reads 
$(PREFIX)_filtered_counts.chrs1-22$(KEEP_CHR).tsv: $(PREFIX)_raw_counts.tsv $(INPUT_MMAP_READS_PILEUP_FILES)
	cat $< | \
	python $(PL)/filter_cnv_sites.py $(PREFIX)_discarded_HetSNVs.tsv $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_rd.tab | \
	python $(PL)/filter_non-autosomal_chr.py $(PREFIX)_discarded_HetSNVs.tsv | \
	python $(PL)/filter_phase_warnings.py $(PREFIX)_discarded_HetSNVs.tsv | \
	python $(PL)/filter_sites_w_mmaps.py $(AMB_MODE) $(PREFIX)_discarded_HetSNVs.tsv $(PREFIX)_sites_w_mmaps.log \
		$(INPUT_MMAP_READS_PILEUP_FILES) > $@

# allelic ratio distrs
$(PREFIX)_ref_allele_ratios.raw_counts.pdf: $(PREFIX)_raw_counts.tsv
	cat $< | python $(PL)/plot_AllelicRatio_distribution.py $(PREFIX)_raw_counts

# counts
$(PREFIX)_raw_counts.tsv: $(INPUT_UNIQ_READS_PILEUP_FILES)
	python $(PL)/pileup2counts.py 1 $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_ref.bed \
	$(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h1.bed $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h2.bed \
	$(PREFIX)_discarded_HetSNVs.tsv \
	$(INPUT_UNIQ_READS_PILEUP_FILES) > $@


