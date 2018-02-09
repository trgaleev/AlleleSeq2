#######################
### PIPELINE PARAMS ###
#######################


### system / executables ##

PL                               := ~/bin/AlleleSeq2
SAMTOOLS                         := ~/bin/samtools-1.3.1/samtools
BEDTOOLS_intersectBed            := ~/bin/bedtools2/bin/intersectBed

### input files / paths ##

REGIONS_FILE                    :=
COUNTS_FILE                       :=
PREFIX                            :=
UNIQ_ALN_FILES                    :=
MMAP_ALN_FILES                    :=

# stats / counts params:

FDR_SIMS                         := 500
FDR_CUTOFF                       := 0.05
Cntthresh_tot                    := 6
Cntthresh_min                    := 0
AMB_MODE                         := adjust # 'adjust' or allelic ratio diff threshold for filtering
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

 
all: $(PREFIX)_region_h1_ratios.filtered_counts.chrs1-22$(KEEP_CHR).pdf $(PREFIX)_region_h1_ratios.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.pdf $(PREFIX)_interesting_regions.FDR-$(FDR_CUTOFF).binom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv $(PREFIX)_interesting_regions.FDR-$(FDR_CUTOFF).betabinom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv

# todo
# this seems to work, but the way it deals with paths, filenames, etc needs to be cleaned up
# here it is even worse
# fix columns
$(PREFIX)_interesting_regions.FDR-$(FDR_CUTOFF).betabinom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv: $(PREFIX)_region_filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv
	awk '{print $$1"\t"$$2"\t"$$3"\t"$$2"\t"$$3"\t0\t0\t"$$2+$$3"\t"$$4"\t1.0"}' $< | \
		sed 's/h1_count\th2_count\t0\t0\t0/cA\tcC\tcG\tcT\tsum_ref_n_alt_cnts/' | \
		sed 's/h1_allele_ratio\t1.0/ref_allele_ratio\tcnv/' > $(PREFIX)_region_filtered_counts_min.$(Cntthresh_tot)-total.$(Cntthresh_min)-lower.tsv.tmp
	Rscript $(PL)/alleledb_calcOverdispersion.R \
		$(PREFIX)_region_filtered_counts_min.$(Cntthresh_tot)-total.$(Cntthresh_min)-lower.tsv.tmp \
		$(PREFIX)_region_FDR-$(FDR_CUTOFF).betabinomial.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min 
	Rscript $(PL)/alleledb_alleleseqBetabinomial.R \
		$(PREFIX)_region_filtered_counts_min.$(Cntthresh_tot)-total.$(Cntthresh_min)-lower.tsv.tmp \
		$(PREFIX)_region_FDR-$(FDR_CUTOFF).betabinomial.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min \
		$(PREFIX)_counts.FDR-$(FDR_CUTOFF).betabinom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv \
		$(PREFIX)_interesting_regions.FDR-$(FDR_CUTOFF).betabinom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv \
		$(PREFIX)_region_FDR-$(FDR_CUTOFF).betabinom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.txt \
		$(FDR_CUTOFF)

# todo fix colnames
$(PREFIX)_interesting_regions.FDR-$(FDR_CUTOFF).binom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv: $(PREFIX)_region_filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv
	python $(PL)/FalsePos_regions.py $< $(FDR_SIMS) $(FDR_CUTOFF) > $(PREFIX)_region_FDR-$(FDR_CUTOFF).binom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.txt
	cat $< | python $(PL)/filter_regions_by_pval.py $(PREFIX)_region_FDR-$(FDR_CUTOFF).binom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.txt > $@

# allelic ratio distrs
$(PREFIX)_region_h1_ratios.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.pdf: $(PREFIX)_region_filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv
	cat $< | python $(PL)/plot_AllelicRatio_distribution.py $(PREFIX)_region_filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt h1_allele_ratio

$(PREFIX)_region_filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv: $(PREFIX)_region_filtered_counts.chrs1-22$(KEEP_CHR).tsv
	cat $< | python $(PL)/filter_regions_by_counts.py $(Cntthresh_tot) $(Cntthresh_min) > $@

# allelic ratio distrs
$(PREFIX)_region_h1_ratios.filtered_counts.chrs1-22$(KEEP_CHR).pdf: $(PREFIX)_region_filtered_counts.chrs1-22$(KEEP_CHR).tsv
	cat $< | python $(PL)/plot_AllelicRatio_distribution.py $(PREFIX)_region_filtered_counts.chrs1-22$(KEEP_CHR) h1_allele_ratio

$(PREFIX)_region_filtered_counts.chrs1-22$(KEEP_CHR).tsv: $(PREFIX)_region_raw_counts_uniq.tsv $(PREFIX)_region_raw_counts_mmap.tsv
	awk '{print $$1"\t"$$2"\t"$$3"\t"$$4"\t"$$5"\t"$$6}' $< | \
	python $(PL)/filter_regions_non-autosomal_chr.py $(REGIONS_FILE) $(PREFIX)_discarded_regions.tsv $(KEEP_CHR) | \
	python $(PL)/filter_regions_w_mmaps.py $(AMB_MODE) \
	$(PREFIX)_discarded_regions.tsv $(PREFIX)_regions_w_mmaps.log $(PREFIX)_region_raw_counts_mmap.tsv > $@

$(PREFIX)_region_raw_counts_mmap.tsv: hets_regions_h2.bed
	cat hets_regions_h1.bed hets_regions_h2.bed | $(BEDTOOLS_intersectBed) -a stdin -b $(MMAP_ALN_FILES) -split -wb -bed | \
	python $(PL)/intersect2counts.py mmap hets_regions_h1.bed > $@

# allowing multiple UNIQ_ALN_FILES files, but probably slower (since became -b instead of -a), also changes to intersect2counts.py to accomodate output file format
#$(PREFIX)_region_counts.tsv: hets_regions_h1.bed hets_regions_h2.bed
#        $(BEDTOOLS_intersectBed) -a $(UNIQ_ALN_FILES) -b hets_regions_h1.bed hets_regions_h2.bed -split -wb -bed | \
#                python $(PL)/intersect2counts.py > $@


# -split -- not counting reads that splice over a het
$(PREFIX)_region_raw_counts_uniq.tsv: hets_regions_h2.bed
	cat hets_regions_h1.bed hets_regions_h2.bed | $(BEDTOOLS_intersectBed) -a stdin -b $(UNIQ_ALN_FILES) -split -wb -bed | \
	python $(PL)/intersect2counts.py uniq hets_regions_h1.bed > $@

hets_regions_h2.bed : $(REGIONS_FILE) $(COUNTS_FILE)
	grep -v '^#chr' $(COUNTS_FILE) | awk '{print $$1"\t"$$2-1"\t"$$2}' | $(BEDTOOLS_intersectBed) -a $(REGIONS_FILE) -b stdin | \
	python $(PL)/refbed2hapcoords.py $(COUNTS_FILE) hets_regions_h1.bed hets_regions_h2.bed

