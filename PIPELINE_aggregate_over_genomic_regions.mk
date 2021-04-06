#######################
### PIPELINE PARAMS ###
#######################


### system / executables ##

PL                               := ~/bin/AlleleSeq2
SAMTOOLS                         := ~/bin/samtools-1.3.1/samtools
BEDTOOLS_intersectBed            := ~/bin/bedtools2/bin/intersectBed

### input files / paths ##

REGIONS_FILE                      :=
COUNTS_FILE                       :=
PREFIX                            :=
UNIQ_ALN_FILES                    :=
MMAP_ALN_FILES                    :=

# stats / counts params:

FDR_SIMS                         := 500
FDR_CUTOFF                       := 0.05
Cntthresh_tot                    := 6
Cntthresh_min                    := 0
#AMB_MODE                         := adjust # 'adjust' or allelic ratio diff threshold for filtering
# only 'adjust' mode only for now: add the 'weaker allele' base counts from multi-mapping reads, if any:
# all of them or until balanced with the stonger to make sure the imbalance is not caused by the multimapping reads 
KEEP_CHR                         := # empty or 'X'


######################
### PIPELINE STEPS ###
######################

empty_string:=

$(info REGIONS_FILE: $(REGIONS_FILE))
$(info COUNTS_FILE: $(COUNTS_FILE))
$(info PREFIX: $(PREFIX))
$(info UNIQ_ALN_FILES: $(UNIQ_ALN_FILES))
$(info MMAP_ALN_FILES: $(MMAP_ALN_FILES))
$(info $(empty_string))



######################
### PIPELINE START ###
######################

 
all:$(PREFIX)_hap1_allele_ratios.region_filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.pdf $(PREFIX)_hap1_allele_ratios.region_filtered_counts.pdf $(PREFIX)_hap1_allele_ratios.region_raw_counts.pdf $(PREFIX)_interesting_regions.FDR-$(FDR_CUTOFF).binom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv $(PREFIX)_interesting_regions.FDR-$(FDR_CUTOFF).betabinom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv


# todo: this seems to work, but the way it deals with paths, filenames, etc needs to be cleaned up
# currently, keeping alleleDB betabinomial scripts with as little modifications as possible
# fix columns and all the workaround with the .tmp files
$(PREFIX)_interesting_regions.FDR-$(FDR_CUTOFF).betabinom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv: $(PREFIX)_region_filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv
	awk '{print $$1"\t"$$2"\t"$$3"\t"$$2"\t"$$3"\t0\t0\t"$$2+$$3"\t"$$4"\t1.0"}' $< | \
		sed 's/hap1_count\thap2_count\t0\t0\t0/cA\tcC\tcG\tcT\tsum_ref_n_alt_cnts/' | \
		sed 's/hap1_allele_ratio\t1.0/ref_allele_ratio\tcnv/' > $(PREFIX)_region_filtered_counts_min.$(Cntthresh_tot)-tot.$(Cntthresh_min)-min.tsv.tmp
	Rscript $(PL)/alleledb_calcOverdispersion.R \
		$(PREFIX)_region_filtered_counts_min.$(Cntthresh_tot)-tot.$(Cntthresh_min)-min.tsv.tmp \
		$(PREFIX)_region_FDR-$(FDR_CUTOFF).betabinomial.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min 
	Rscript $(PL)/alleledb_alleleseqBetabinomial.R \
		$(PREFIX)_region_filtered_counts_min.$(Cntthresh_tot)-tot.$(Cntthresh_min)-min.tsv.tmp \
		$(PREFIX)_region_FDR-$(FDR_CUTOFF).betabinomial.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min \
		$(PREFIX)_counts.FDR-$(FDR_CUTOFF).betabinom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv \
		$(PREFIX)_interesting_regions.FDR-$(FDR_CUTOFF).betabinom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tmp \
		$(PREFIX)_region_FDR-$(FDR_CUTOFF).betabinom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.txt \
		$(FDR_CUTOFF)
	awk '{print $$1"\t"$$2"\t"$$3"\t"$$8"\t"$$9"\t"$$11}' \
		$(PREFIX)_interesting_regions.FDR-$(FDR_CUTOFF).betabinom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tmp | \
		sed 's/sum_ref_n_alt_cnts/hap1_allele_ratio/' > $@ 
	#rm $(PREFIX)_region_filtered_counts_min.$(Cntthresh_tot)-tot.$(Cntthresh_min)-min.tsv.tmp \
	#	$(PREFIX)_region_FDR-$(FDR_CUTOFF).betabinom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tmp

$(PREFIX)_interesting_regions.FDR-$(FDR_CUTOFF).binom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv: $(PREFIX)_region_filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv
	python $(PL)/FalsePos_regions.py $< $(FDR_SIMS) $(FDR_CUTOFF) > $(PREFIX)_region_FDR-$(FDR_CUTOFF).binom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.txt
	cat $< | python $(PL)/filter_regions_by_pval.py $(PREFIX)_region_FDR-$(FDR_CUTOFF).binom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.txt > $@

# allelic ratio distrs
$(PREFIX)_hap1_allele_ratios.region_filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.pdf: $(PREFIX)_region_filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv
	Rscript $(PL)/plot_AllelicRatio_distribution.R $< $(PREFIX) region_filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt hap1_allele_ratio

# filter based on total counts and min per hap count
# and in non-autosomal chr, optionally keeping X;
$(PREFIX)_region_filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv: $(PREFIX)_region_filtered_counts.tsv
	cat $< | \
	python $(PL)/filter_regions_non-autosomal_chr.py $(REGIONS_FILE) $(KEEP_CHR) | \
	python $(PL)/filter_regions_by_counts.py $(Cntthresh_tot) $(Cntthresh_min) > $@

# allelic ratio distrs
$(PREFIX)_hap1_allele_ratios.region_filtered_counts.pdf: $(PREFIX)_region_filtered_counts.tsv
	Rscript $(PL)/plot_AllelicRatio_distribution.R $< $(PREFIX) region_filtered_counts hap1_allele_ratio

$(PREFIX)_region_filtered_counts.tsv: $(PREFIX)_region_raw_counts_uniq.tsv $(PREFIX)_region_raw_counts_mmap.tsv
	cat $< | \
#	python $(PL)/filter_regions_w_mmaps.py $(AMB_MODE) \
	python $(PL)/filter_regions_w_mmaps.py adjust \
	$(PREFIX)_discarded_regions_warn-mmaps.log $(PREFIX)_mmap_reads_over_regions.log $(PREFIX)_region_raw_counts_mmap.tsv > $@

# this can be memory-intensive and fail without error message (empty $(PREFIX)_region_raw_counts_mmap.tsv file)
$(PREFIX)_region_raw_counts_mmap.tsv: hets_regions_hap2.bed
	cat hets_regions_hap1.bed hets_regions_hap2.bed | $(BEDTOOLS_intersectBed) -a stdin -b $(MMAP_ALN_FILES) -split -wb -bed | \
	python $(PL)/intersect2counts.py mmap hets_regions_hap1.bed > $@

# allelic ratio distrs
$(PREFIX)_hap1_allele_ratios.region_raw_counts.pdf: $(PREFIX)_region_raw_counts_uniq.tsv
	Rscript $(PL)/plot_AllelicRatio_distribution.R $< $(PREFIX) region_raw_counts hap1_allele_ratio



# -split -- not counting reads that splice over a het
$(PREFIX)_region_raw_counts_uniq.tsv: hets_regions_hap2.bed
	cat hets_regions_hap1.bed hets_regions_hap2.bed | $(BEDTOOLS_intersectBed) -a stdin -b $(UNIQ_ALN_FILES) -split -wb -bed | \
	python $(PL)/intersect2counts.py uniq hets_regions_hap1.bed > $@

hets_regions_hap2.bed : $(REGIONS_FILE) $(COUNTS_FILE)
	grep -v '^#chr' $(COUNTS_FILE) | awk '{print $$1"\t"$$2-1"\t"$$2}' | $(BEDTOOLS_intersectBed) -a $(REGIONS_FILE) -b stdin | \
	python $(PL)/refbed2hapcoords.py $(COUNTS_FILE) hets_regions_hap1.bed hets_regions_hap2.bed

