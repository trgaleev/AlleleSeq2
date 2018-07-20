USAGE := "make makePersonalGenome VCF_SAMPLE_ID=sample_name FILE_PATH_VCF=path_to_vcf FILE_PATH_BAM=path_to_wgs_bam"

#### executables, system params
VCF2DIPLOID_DIR               := ~/bin/vcf2diploid_v0.2.6a/
JAVA                          := /usr/bin/java
BOWTIE_build                  := ~/bin/bowtie-1.1.1/bowtie-build
BOWTIE2_build                 := ~/bin/bowtie2-2.3.1/bowtie2-build
LIFTOVER                      := ~/bin/liftOver
BEDTOOLS_intersectBed         := ~/bin/bedtools2/bin/intersectBed
BEDTOOLS_coverageBed          := ~/bin/bedtools2/bin/coverageBed
SAMTOOLS                      := ~/bin/samtools-1.3.1/samtools
MAX_RAM                       := 45G
PL                            := ~/bin/AlleleSeq2/

# matters if used STAR genomeGenerate only
STAR                          := ~/bin/STAR/bin/Linux_x86_64/STAR
N_THREADS                     := 1        
STAR_limitGenomeGenerateRAM   := 31000000000 # (bytes, default in STAR)

# aligner:
#  'bowtie1' or 'STAR'; 

ALIGNER                       := STAR

#### paths

REFGENOME_VERSION             := GRCh37   #GRCh38 or CRCh37

ifeq ($(REFGENOME_VERSION), GRCh38)

## Primary assembly

  ### Assebmled chromosomes (only 1-22, X, Y, M will be used by vcf2diploid if the reference has other contigs) with PAR and repeat array regions hard masked:
  REFGENOME        := ~/refs_annotations/UCSC_Analysis_Set_GRCh38/hg38.analysisSet.fasta
  ### Unlocalized sequences (_random) and unplaced sequences (chrU_) of the primary assebmly (in addition to assembled chromosomes)
  ### + EBV&decoy contigs (without alternate contigs, alternate scaffolds or alternate loci (_alt)).
  ### These will be added to the diploid personal genome to siphon off reads corresponding to those during mapping
  ### will not be used for readcounting, AS calling, etc
  ADDNL_SEQNS      := ~/refs_annotations/UCSC_Analysis_Set_GRCh38/hg38.analysisSet_unplaced_unlocalized_EBV_contigs_only.fasta

  ANNOTATION       := ~/refs_annotations/gencode.v24.annotation.gtf

else ifeq ($(REFGENOME_VERSION), GRCh37)

  REFGENOME        := ~/refs_annotations/hg19_ucsc.fasta
  ADDNL_SEQNS      := ~/refs_annotations/hg19_ucsc_non_chr_scaffolds_only.fasta
  ANNOTATION       := ~/refs_annotations/gencode.v19.annotation.gtf

endif

ALLOSOMES_M        := X,Y
ALLOSOMES_P        := X,Y
# if globally phased and the vcf follows GT=pat|mat convention, can be specified accordingly for male/female
# otherwise will use all combinations to make sure all possible seqs
# are there when mapping; hets from non-autosomal should be removed from AS calls then




# (required) inputs ####
VCF_SAMPLE_ID   := NULL
FILE_PATH_BAM   := NULL
FILE_PATH_VCF   := NULL  # if separate vcfs, use this for SNVs and the following for indels and svs 

# specify only if using separate vcf files
FILE_PATH_VCF_INDELS := 
FILE_PATH_VCF_SVS    :=   


OUTPUT_DIR := pgenome_$(VCF_SAMPLE_ID)


#### aligner idx
empty_string  :=
ifeq ($(ALIGNER),bowtie1)
  PGENOME_IDX_TARGET = $(OUTPUT_DIR)/bowtie_build_diploid.log
else ifeq ($(ALIGNER),STAR)
  PGENOME_IDX_TARGET = $(OUTPUT_DIR)/STAR_idx_diploid_Log.out
else ifeq ($(ALIGNER)$(ADDNL_SEQNS),bowtie2$(empty_string))
  PGENOME_IDX_TARGET = $(OUTPUT_DIR)/bowtie2_build_diploid.log
  BOWTIE2_SEQ_IN = $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hap1.fa,$(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hap2.fa 
  $(info BOWTIE2_SEQ_IN: $(BOWTIE2_SEQ_IN))
else ifeq ($(ALIGNER),bowtie2)
  PGENOME_IDX_TARGET = $(OUTPUT_DIR)/bowtie2_build_diploid.log
  BOWTIE2_SEQ_IN = $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hap1.fa,$(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hap2.fa,$(ADDNL_SEQNS)
  $(info BOWTIE2_SEQ_IN: $(BOWTIE2_SEQ_IN))
endif


$(info REFGENOME_VERSION: $(REFGENOME_VERSION))
$(info REFGENOME: $(REFGENOME))
$(info ANNOTATION: $(ANNOTATION))
$(info ADDNL_SEQNS: $(ADDNL_SEQNS))
$(info FILE_PATH_VCF: $(FILE_PATH_VCF))
$(info FILE_PATH_VCF_INDELS: $(FILE_PATH_VCF_INDELS))
$(info FILE_PATH_VCF_SVS: $(FILE_PATH_VCF_SVS))
$(info FILE_PATH_BAM: $(FILE_PATH_BAM))
$(info OUTPUT_DIR: $(OUTPUT_DIR))
$(info PGENOME_IDX_TARGET: $(PGENOME_IDX_TARGET))


#### pipeline start ######

all: $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_hap1.bed $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_hap2.bed $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_rd.tab $(PGENOME_IDX_TARGET) 


# vcf2diploid

$(OUTPUT_DIR)/maternal.chain $(OUTPUT_DIR)/paternal.chain: $(FILE_PATH_VCF)
	@echo -e "$(USAGE)"
	mkdir -p $(OUTPUT_DIR)
	$(JAVA) -Xmx$(MAX_RAM) -jar $(VCF2DIPLOID_DIR)/vcf2diploid.jar \
	-id $(VCF_SAMPLE_ID) \
	-pass \
	-chr $(REFGENOME) \
	-vcf $(FILE_PATH_VCF_SVS) $(FILE_PATH_VCF_INDELS) $(FILE_PATH_VCF) \
	-outDir $(OUTPUT_DIR)



# prepare .bed files w all hetSNVs in both coord systems
# if chr name convention is different in .vcf and $(REFGENOME) use the one from $(REFGENOME) 
# also make it more general: substitute 'paternal' with 'hap1' and 'maternal' with 'hap2' to match GT, hap1|hap2, in the vcf

$(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_ref.bed: $(FILE_PATH_VCF) $(OUTPUT_DIR)/maternal.chain $(OUTPUT_DIR)/paternal.chain
	cat $(FILE_PATH_VCF) | python $(PL)/get_hetSNV_bed.py $(VCF_SAMPLE_ID) $(REFGENOME) > $@

$(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_hap2.bed: $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_ref.bed $(OUTPUT_DIR)/maternal.chain
	sed 's/maternal/hap2/g' $(OUTPUT_DIR)/maternal.chain | $(LIFTOVER) $< stdin $@ $(subst .bed,,$@).not_lifted.bed

$(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_hap1.bed: $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_ref.bed $(OUTPUT_DIR)/paternal.chain
	sed 's/paternal/hap1/g' $(OUTPUT_DIR)/paternal.chain | $(LIFTOVER) $< stdin $@ $(subst .bed,,$@).not_lifted.bed



#$(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hap2.fa: $(OUTPUT_DIR)/maternal.chain $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_ref.bed
#	cat $(shell awk  '{print $$1"_maternal.fa"}' $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_ref.bed | uniq) | \
#	sed 's/maternal/hap2/g' > $@

#$(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hap1.fa: $(OUTPUT_DIR)/paternal.chain $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_ref.bed
#	cat $(shell awk  '{print $$1"_paternal.fa"}' $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_ref.bed | uniq) | \
#	sed 's/paternal/hap1/g' > $@
	

$(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hap2.fa: $(OUTPUT_DIR)/maternal.chain $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_ref.bed
	cat $(shell ls $(OUTPUT_DIR)/*[0-9]_$(VCF_SAMPLE_ID)_maternal.fa | sort -V)  \
	    $(shell ls $(OUTPUT_DIR)/*[$(ALLOSOMES_M)]_$(VCF_SAMPLE_ID)_maternal.fa)  \
	    $(shell ls $(OUTPUT_DIR)/*[M,MT]_$(VCF_SAMPLE_ID)_maternal.fa) | \
	sed 's/maternal/hap2/g' > $@

$(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hap1.fa: $(OUTPUT_DIR)/paternal.chain $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_ref.bed
	cat $(shell ls $(OUTPUT_DIR)/*[0-9]_$(VCF_SAMPLE_ID)_paternal.fa | sort -V)  \
	    $(shell ls $(OUTPUT_DIR)/*[$(ALLOSOMES_P)]_$(VCF_SAMPLE_ID)_paternal.fa) | \
	sed 's/paternal/hap1/g' > $@

# if star
# will generate lifted annotation files, but will not use them when generating STAR genome indices, since
# star v 2.4.1a can do this on the fly; thus --sjdbGTFfile and --sjdbOverhang should be specified when mapping the reads

$(OUTPUT_DIR)/maternal.$(notdir $(ANNOTATION)): $(OUTPUT_DIR)/maternal.chain 
	cat $(ANNOTATION) | $(LIFTOVER) -gff stdin $< $@ $(OUTPUT_DIR)/maternal.notLifted

$(OUTPUT_DIR)/paternal.$(notdir $(ANNOTATION)): $(OUTPUT_DIR)/paternal.chain 
	cat $(ANNOTATION) | $(LIFTOVER) -gff stdin $< $@ $(OUTPUT_DIR)/paternal.notLifted

$(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_diploid.$(notdir $(ANNOTATION)): $(OUTPUT_DIR)/maternal.$(notdir $(ANNOTATION)) $(OUTPUT_DIR)/paternal.$(notdir $(ANNOTATION))
	cat $(OUTPUT_DIR)/maternal.$(notdir $(ANNOTATION)) $(OUTPUT_DIR)/paternal.$(notdir $(ANNOTATION)) | sed 's/maternal/hap2/g' | sed 's/paternal/hap1/g'  > $@


# star idx

$(OUTPUT_DIR)/STAR_idx_diploid_Log.out: $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_diploid.$(notdir $(ANNOTATION)) $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hap2.fa $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hap1.fa
	mkdir -p $(OUTPUT_DIR)/STAR_idx_diploid
	$(STAR) \
		--runThreadN $(N_THREADS) \
		--limitGenomeGenerateRAM $(STAR_limitGenomeGenerateRAM) \
		--runMode genomeGenerate \
		--genomeDir $(OUTPUT_DIR)/STAR_idx_diploid \
		--genomeFastaFiles $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hap1.fa $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hap2.fa $(ADDNL_SEQNS) \
		--outFileNamePrefix $(OUTPUT_DIR)/STAR_idx_diploid_


### attempting to make this faster and deal easier with differences in chr naming conventions:

#$(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp.bed.counts: $(OUTPUT_DIR)/maternal.chain $(OUTPUT_DIR)/paternal.chain
#	$(VCF2DIPLOID_DIR)/vcf2snp -c $(VCF_SAMPLE_ID) -p 1 -r 1 -s 1 $(FILE_PATH_VCF) > $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp
#	grep -w '^1\|^2\|^3\|^4\|^5\|^6\|^7\|^8\|^9\|^10\|^11\|^12\|^13\|^14\|^15\|^16\|^17\|^18\|^19\|^20\|^21\|^22\|^X\|^Y\|^MT' $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp | \
#		awk '{print $$1"\t"($$2-1000)"\t"($$2+1000)"\t"$$1"_"($$2-1000)"_"($$2+1000)}' | grep -v "-" > $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp.bed
#	$(BEDTOOLS_intersectBed) -wo -bed -a $(FILE_PATH_BAM) -b $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp.bed | awk '{print $$16}' | sort | uniq -c > $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp.bed.counts


#$(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp.bed.mediancount: $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp.bed.counts
#	cat $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp.bed.counts | sort -n -k 1 | awk '{ lines[NR]=$$0; } END { print lines[int(NR/2)+1] }' | awk '{print $$1}' > $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp.bed.mediancount


#$(OUTPUT_DIR)/$(VCF_SAMPLE_ID).alleleSeqInput.snp $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).alleleSeqInput.cnv: $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp.bed.mediancount
#	cat $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp.bed.counts | awk '{getline avg<"$(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp.bed.mediancount"; print $$2"\t"$$1/avg }' > $(OUTPUT_DIR)/tmp.cnv
#	cat $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp | awk '{print $$1"_"($$2-1000)"_"($$2+1000)"\t"$$0}' | sort -k 1 > $(OUTPUT_DIR)/tmp.snp
#	awk 'NR==FNR {h[$$1] = $$0; next} {print h[$$1]"\t"$$0}' $(OUTPUT_DIR)/tmp.snp $(OUTPUT_DIR)/tmp.cnv > $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snpANDcnv
#	cat $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snpANDcnv | awk '{print $$2"\t"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7"\t"$$8}' > $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).alleleSeqInput.snp
#	echo -e "chrm\tsnppos\trd" > $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).alleleSeqInput.cnv
#	cat $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snpANDcnv | awk '{print $$2"\t"$$3"\t"$$10}' >> $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).alleleSeqInput.cnv
#	rm $(OUTPUT_DIR)/tmp.cnv $(OUTPUT_DIR)/tmp.snp



# this seems to reproduce previous approach, but ~3 times faster
# requires chr in the .bed to be in the same order as in the .bam file
# #todo: assess if should be using (by prefiltering with samtools) only non-duplicated reads
# in the wgs bam file
# also test samtools bedcov

#$(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_rd.tab: $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_ref.bed
#	$(SAMTOOLS) view -H $(FILE_PATH_BAM) | grep -P "@SQ\tSN:" | awk -F"\t" '{print $$1"\t"$$2"\t"$$3}' | \
#	sed 's/@SQ\tSN://' | sed 's/\tLN:/\t/' > $(OUTPUT_DIR)/genome.txt 
#	awk '{print $$1"\t"($$3-1000)"\t"($$3+1000)"\t"$$4}' $< | grep -v "-" | \
#	python $(PL)/feed_sorted_check_names.py $(OUTPUT_DIR)/genome.txt | \
#	$(BEDTOOLS_coverageBed) -a stdin -b $(FILE_PATH_BAM) \
#		-g $(OUTPUT_DIR)/genome.txt -sorted -counts | python $(PL)/calculate_rd.py > $@
	


# actually, having the larger file as -a is better for the memory

$(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_rd.tab: $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_ref.bed
	$(SAMTOOLS) view -H $(FILE_PATH_BAM) | grep -P "@SQ\tSN:" | awk -F"\t" '{print $$1"\t"$$2"\t"$$3}' | \
	sed 's/@SQ\tSN://' | sed 's/\tLN:/\t/' > $(OUTPUT_DIR)/genome.txt 
	awk '{print $$1"\t"($$3-1000)"\t"($$3+1000)"\t"$$4}' $< | grep -v "-" | \
	python $(PL)/feed_sorted_check_names.py $(OUTPUT_DIR)/genome.txt | \
	$(BEDTOOLS_intersectBed) -a $(FILE_PATH_BAM) -b stdin \
		-g $(OUTPUT_DIR)/genome.txt -sorted -wb -bed | python $(PL)/calculate_rd.py > $@


#$(OUTPUT_DIR)/bowtie_build_diploid.log: $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hap2.fa $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hap1.fa 
#	$(BOWTIE_build) --offrate 2 \
#		$(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hap2.fa,$(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hap1.fa $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_diploid > $@



$(OUTPUT_DIR)/bowtie2_build_diploid.log: $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hap1.fa $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hap2.fa 
	$(BOWTIE2_build) --threads $(N_THREADS) $(BOWTIE2_SEQ_IN) $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_diploid > $@
