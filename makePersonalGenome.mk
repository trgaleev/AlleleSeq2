USAGE := "make makePersonalGenome VCF_SAMPLE_ID=sample_name FILE_PATH_VCF=path_to_vcf FILE_PATH_BAM=path_to_wgs_bam"

#### executables, system params
VCF2DIPLOID_DIR               := ~/bin/vcf2diploid_v0.2.6a/
JAVA                          := /usr/bin/java
#BOWTIE_build                  := ~/bin/bowtie-1.1.1/bowtie-build
#BOWTIE2_build                 := ~/bin/bowtie2-2.3.1/bowtie2-build
LIFTOVER                      := ~/bin/liftOver
BEDTOOLS_intersectBed         := ~/bin/bedtools2/bin/intersectBed
SAMTOOLS                      := ~/bin/samtools-1.3.1/samtools
MAX_RAM                       := 45G
PL                            := ~/bin/AlleleSeq2/

# matters if used STAR genomeGenerate only
STAR                          := ~/bin/STAR/STAR-2.6.0c/bin/Linux_x86_64/STAR
N_THREADS                     := 1        
STAR_limitGenomeGenerateRAM   := 31000000000 # (bytes, default in STAR)

# aligner:
#  'bowtie1' or 'STAR'; Only using STAR now 

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



$(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_rd.tab: $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_ref.bed
	$(SAMTOOLS) view -H $(FILE_PATH_BAM) | grep -P "@SQ\tSN:" | awk -F"\t" '{print $$1"\t"$$2"\t"$$3}' | \
	sed 's/@SQ\tSN://' | sed 's/\tLN:/\t/' > $(OUTPUT_DIR)/genome.txt 
	awk '{print $$1"\t"($$3-1000)"\t"($$3+1000)"\t"$$4}' $< | grep -v "-" | \
	python $(PL)/feed_sorted_check_names.py $(OUTPUT_DIR)/genome.txt | \
	$(BEDTOOLS_intersectBed) -a $(FILE_PATH_BAM) -b stdin \
		-g $(OUTPUT_DIR)/genome.txt -sorted -wb -bed | python $(PL)/calculate_rd.py > $@


#$(OUTPUT_DIR)/bowtie2_build_diploid.log: $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hap1.fa $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hap2.fa 
#	$(BOWTIE2_build) --threads $(N_THREADS) $(BOWTIE2_SEQ_IN) $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_diploid > $@
