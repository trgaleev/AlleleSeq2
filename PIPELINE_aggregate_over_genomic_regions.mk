#######################
### PIPELINE PARAMS ###
#######################


### system / executables ##

PL                               := ~/bin/AlleleSeq2
SAMTOOLS                         := ~/bin/samtools-1.3.1/samtools
BEDTOOLS_intersectBed            := ~/bin/bedtools2/bin/intersectBed

### input files / paths ##

INTERVALS_FILE                    :=
COUNTS_FILE                       :=
UNIQ_ALN_FILE                     :=

# stats / counts params:

FDR_SIMS                         := 500
FDR_CUTOFF                       := 0.05

######################
### PIPELINE STEPS ###
######################

empty_string:=

$(info INTERVALS_FILE: $(INTERVALS_FILE))
$(info COUNTS_FILE: $(COUNTS_FILE))
$(info $(empty_string))



######################
### PIPELINE START ###
######################

 
all: interval_interestingHets.binom.tsv interval_interestingHets.betabinom.tsv

# this seems to work, but the way it deals with paths, filenames, etc needs to be cleaned up
# here it is even worse
interval_interestingHets.betabinom.tsv: interval_counts.tsv
	awk '{print $$1"\t"$$2"\t"$$3"\t"$$2"\t"$$3"\t0\t0\t"$$2+$$3"\t"$$4"\t1.0"}' interval_counts.tsv | sed 's/h1_count\th2_count\t0\t0\t0/cA\tcC\tcG\tcT\tsum_ref_n_alt_cnts/' | sed 's/h1_allele_ratio\t1.0/ref_allele_ratio\tcnv/' > interval_counts.tsv.tmp 
	Rscript $(PL)/alleledb_calcOverdispersion.R interval_counts.tsv.tmp interval_betabinomial 
	Rscript $(PL)/alleledb_alleleseqBetabinomial.R interval_counts.tsv.tmp interval_betabinomial interval_counts.betabinom.tsv interval_interestingHets.betabinom.tsv interval_FDR.betabinomial.txt $(FDR_CUTOFF)

interval_interestingHets.binom.tsv: interval_counts.tsv
	python $(PL)/FalsePos_interval.py $< $(FDR_SIMS) $(FDR_CUTOFF) > interval_FDR.binom.txt
	cat $< | python $(PL)/filter_by_pval_interval.py interval_FDR.binom.txt > $@

interval_counts_h1_allele_ratios.pdf: interval_counts.tsv
	cat $< | python $(PL)/plot_AllelicRatio_distribution.py $< h1_allele_ratio

interval_counts.tsv: hets_intervals_h1.bed hets_intervals_h2.bed
	$(BEDTOOLS_intersectBed) -a $(UNIQ_ALN_FILE) -b hets_intervals_h1.bed hets_intervals_h2.bed -split -wb -bed | \
		python $(PL)/intersect2counts.py > $@

hets_intervals_h1.bed hets_intervals_h2.bed: $(INTERVALS_FILE) $(COUNTS_FILE)
	grep -v '^#chr' $(COUNTS_FILE) | awk '{print $$1"\t"$$2-1"\t"$$2}' | $(BEDTOOLS_intersectBed) -a $(INTERVALS_FILE) -b stdin | \
		python $(PL)/refbed2hapcoords.py $(COUNTS_FILE) hets_intervals_h1.bed hets_intervals_h2.bed

