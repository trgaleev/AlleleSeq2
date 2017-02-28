USAGE := "make makePersonalGenome VCF_SAMPLE_ID=sample_name FILE_PATH_VCF=path_to_vcf FILE_PATH_BAM=path_to_wgs_bam"


#### executables, system params
VCF2DIPLOID_DIR               := ~/bin/vcf2diploid_v0.2.6a/
JAVA                          := /usr/bin/java
BOWTIE_build                  := ~/bin/bowtie-1.1.1/bowtie-build
LIFTOVER                      := ~/bin/liftOver
BEDTOOLS_intersectBed         := ~/bin/bedtools2/bin/intersectBed
MAX_RAM                       := 45G
PL                            := ~/bin/pgenome_mkfiles/

# matters if used STAR genomeGenerate only
STAR                          := ~/bin/STAR/bin/Linux_x86_64/STAR
N_THREADS                     := 1        

# aligner:
#  'bowtie1' as in AlleleSeq v1.2a and before or 
#  'STAR' for spliced alignment of RNA-seq data 
# and type of genome idx:
#  'separate' maternal and paternal genome indices separately (AlleleSeq v1.2a)
#  'combined' one combined diploid genome idx with all maternal and paternal chromosomes 

ALIGNER                       := STAR
GENOME_IDX_TYPE               := combined

#### paths

# hg38
#REFGENOME        := ~/refs_annotations/Homo_sapiens_assembly38.fasta
#ANNOTATION       := ~/refs_annotations/gencode.v24.annotation.gtf
# hg19
REFGENOME        := ~/refs_annotations/human_g1k_v37_decoy.fasta
ANNOTATION       := ~/refs_annotations/gencode.v19.annotation.gtf


# (required) inputs ####
FILE_PATH_BAM   := NULL
FILE_PATH_VCF   := NULL  # if separate, use this for SNVs and the following for indels and svs 
VCF_SAMPLE_ID   := NULL

# specify only if using separate vcf files
FILE_PATH_VCF_INDELS := 
FILE_PATH_VCF_SVS    :=   


OUTPUT_DIR := pgenome_$(VCF_SAMPLE_ID)


#### type of genome idx

ifeq ($(GENOME_IDX_TYPE),separate)
 SUFFIX = paternal
else ifeq ($(GENOME_IDX_TYPE),combined)
 SUFFIX = diploid
endif

#### bowtie vs STAR

ifeq ($(ALIGNER),bowtie1)
  GENOME_IDX_TARGET = $(OUTPUT_DIR)/bowtie_build_$(SUFFIX).log
else ifeq ($(ALIGNER),STAR)
  GENOME_IDX_TARGET = $(OUTPUT_DIR)/STAR_idx_$(SUFFIX)_Log.out
endif





#$(warning GENOME_IDX_TARGET: $(GENOME_IDX_TARGET))




#### pipeline start ######

all: $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).alleleSeqInput.snp $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).alleleSeqInput.cnv $(GENOME_IDX_TARGET)  




$(OUTPUT_DIR)/maternal.chain $(OUTPUT_DIR)/paternal.chain: $(FILE_PATH_VCF)  
	@echo -e "$(USAGE)"
	mkdir $(OUTPUT_DIR)
	$(JAVA) -Xmx$(MAX_RAM) -jar $(VCF2DIPLOID_DIR)/vcf2diploid.jar -id $(VCF_SAMPLE_ID) -pass -chr $(REFGENOME) -vcf $(FILE_PATH_VCF_SVS) $(FILE_PATH_VCF_INDELS) $(FILE_PATH_VCF) -outDir $(OUTPUT_DIR)



$(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp.bed.counts: $(OUTPUT_DIR)/maternal.chain $(OUTPUT_DIR)/paternal.chain
	$(VCF2DIPLOID_DIR)/vcf2snp -c $(VCF_SAMPLE_ID) -p 1 -r 1 -s 1 $(FILE_PATH_VCF) > $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp
	grep -w '^1\|^2\|^3\|^4\|^5\|^6\|^7\|^8\|^9\|^10\|^11\|^12\|^13\|^14\|^15\|^16\|^17\|^18\|^19\|^20\|^21\|^22\|^X\|^Y\|^MT' $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp | \
		awk '{print $$1"\t"($$2-1000)"\t"($$2+1000)"\t"$$1"_"($$2-1000)"_"($$2+1000)}' | grep -v "-" > $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp.bed
	$(BEDTOOLS_intersectBed) -wo -bed -a $(FILE_PATH_BAM) -b $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp.bed | awk '{print $$16}' | sort | uniq -c > $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp.bed.counts


$(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp.bed.mediancount: $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp.bed.counts
	cat $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp.bed.counts | sort -n -k 1 | awk '{ lines[NR]=$$0; } END { print lines[int(NR/2)+1] }' | awk '{print $$1}' > $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp.bed.mediancount


$(OUTPUT_DIR)/$(VCF_SAMPLE_ID).alleleSeqInput.snp $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).alleleSeqInput.cnv: $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp.bed.mediancount
	cat $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp.bed.counts | awk '{getline avg<"$(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp.bed.mediancount"; print $$2"\t"$$1/avg }' > $(OUTPUT_DIR)/tmp.cnv
	cat $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snp | awk '{print $$1"_"($$2-1000)"_"($$2+1000)"\t"$$0}' | sort -k 1 > $(OUTPUT_DIR)/tmp.snp
	awk 'NR==FNR {h[$$1] = $$0; next} {print h[$$1]"\t"$$0}' $(OUTPUT_DIR)/tmp.snp $(OUTPUT_DIR)/tmp.cnv > $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snpANDcnv
	cat $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snpANDcnv | awk '{print $$2"\t"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7"\t"$$8}' > $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).alleleSeqInput.snp
	echo -e "chrm\tsnppos\trd" > $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).alleleSeqInput.cnv
	cat $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).snpANDcnv | awk '{print $$2"\t"$$3"\t"$$10}' >> $(OUTPUT_DIR)/$(VCF_SAMPLE_ID).alleleSeqInput.cnv
	rm $(OUTPUT_DIR)/tmp.cnv $(OUTPUT_DIR)/tmp.snp


$(OUTPUT_DIR)/bowtie_build_maternal.log: $(OUTPUT_DIR)/maternal.chain $(OUTPUT_DIR)/paternal.chain
	$(BOWTIE_build) --offrate 2 \
		$(OUTPUT_DIR)/1_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/2_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/3_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/4_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/5_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/6_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/7_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/8_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/9_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/10_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/11_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/12_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/13_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/14_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/15_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/16_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/17_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/18_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/19_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/20_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/21_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/22_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/X_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/Y_$(VCF_SAMPLE_ID)_maternal.fa,$(OUTPUT_DIR)/MT_$(VCF_SAMPLE_ID)_maternal.fa $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_maternal > $@


$(OUTPUT_DIR)/bowtie_build_paternal.log: $(OUTPUT_DIR)/bowtie_build_maternal.log
	$(BOWTIE_build) --offrate 2 \
		$(OUTPUT_DIR)/1_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/2_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/3_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/4_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/5_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/6_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/7_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/8_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/9_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/10_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/11_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/12_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/13_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/14_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/15_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/16_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/17_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/18_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/19_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/20_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/21_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/22_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/X_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/Y_$(VCF_SAMPLE_ID)_paternal.fa,$(OUTPUT_DIR)/MT_$(VCF_SAMPLE_ID)_paternal.fa $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_paternal > $@




# if star
# will generate lifted annotation files, but will not use them when generating STAR genome indices, since
# star v 2.4.1a can do this on the fly, thus --sjdbGTFfile and --sjdbOverhang should be specified when mapping the reads


$(OUTPUT_DIR)/maternal.$(notdir $(ANNOTATION)): $(OUTPUT_DIR)/maternal.chain 
	sed 's/^chr//g' $(ANNOTATION) | $(LIFTOVER) -gff stdin $< $@ $(OUTPUT_DIR)/maternal.notLifted

$(OUTPUT_DIR)/paternal.$(notdir $(ANNOTATION)): $(OUTPUT_DIR)/paternal.chain 
	sed 's/^chr//g' $(ANNOTATION) | $(LIFTOVER) -gff stdin $< $@ $(OUTPUT_DIR)/paternal.notLifted


# for combined, will also make it more general: substitute 'paternal' with h1 and 'maternal' with h2 to match GT, h1|h2, in the vcf

$(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_diploid.$(notdir $(ANNOTATION)): $(OUTPUT_DIR)/maternal.$(notdir $(ANNOTATION)) $(OUTPUT_DIR)/paternal.$(notdir $(ANNOTATION))
	cat $(OUTPUT_DIR)/maternal.$(notdir $(ANNOTATION)) $(OUTPUT_DIR)/paternal.$(notdir $(ANNOTATION)) | sed 's/maternal/h2/g' | sed 's/paternal/h1/g'  > $@


# also prepare .bed files w all hetSNVs in both coord systems, will be useful for combined

$(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_ref.bed: $(FILE_PATH_VCF)
	cat $(FILE_PATH_VCF) | python $(PL)/get_hetSNV_bed.py $(VCF_SAMPLE_ID) > $@


$(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h2.bed: $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_ref.bed
	sed 's/maternal/h2/g' $(OUTPUT_DIR)/maternal.chain | $(LIFTOVER) $< stdin $@ $(subst .bed,,$@).not_lifted.bed

$(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h1.bed: $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_ref.bed
	sed 's/paternal/h1/g' $(OUTPUT_DIR)/paternal.chain | $(LIFTOVER) $< stdin $@ $(subst .bed,,$@).not_lifted.bed



$(OUTPUT_DIR)/STAR_idx_maternal_Log.out: $(OUTPUT_DIR)/maternal.$(notdir $(ANNOTATION))
	mkdir $(OUTPUT_DIR)/STAR_idx_maternal
	$(STAR) \
		--runThreadN $(N_THREADS) \
		--runMode genomeGenerate \
		--genomeDir $(OUTPUT_DIR)/STAR_idx_maternal \
		--genomeFastaFiles \
			$(OUTPUT_DIR)/1_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/2_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/3_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/4_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/5_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/6_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/7_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/8_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/9_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/10_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/11_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/12_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/13_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/14_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/15_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/16_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/17_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/18_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/19_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/20_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/21_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/22_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/X_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/Y_$(VCF_SAMPLE_ID)_maternal.fa \
			$(OUTPUT_DIR)/MT_$(VCF_SAMPLE_ID)_maternal.fa \
		--outFileNamePrefix $(OUTPUT_DIR)/STAR_idx_maternal_


$(OUTPUT_DIR)/STAR_idx_paternal_Log.out: $(OUTPUT_DIR)/paternal.$(notdir $(ANNOTATION)) $(OUTPUT_DIR)/STAR_idx_maternal_Log.out
	mkdir $(OUTPUT_DIR)/STAR_idx_paternal
	$(STAR) \
		--runThreadN $(N_THREADS) \
		--runMode genomeGenerate \
		--genomeDir $(OUTPUT_DIR)/STAR_idx_paternal \
		--genomeFastaFiles \
			$(OUTPUT_DIR)/1_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/2_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/3_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/4_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/5_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/6_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/7_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/8_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/9_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/10_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/11_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/12_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/13_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/14_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/15_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/16_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/17_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/18_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/19_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/20_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/21_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/22_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/X_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/Y_$(VCF_SAMPLE_ID)_paternal.fa \
			$(OUTPUT_DIR)/MT_$(VCF_SAMPLE_ID)_paternal.fa \
		--outFileNamePrefix $(OUTPUT_DIR)/STAR_idx_paternal_



# if combined, will use chrs 1-22 as well as X,Y, MT from both haplotypes as single diploid genome, even though it doesn't make much sense to use all of the latter. 
# such as paternal MT, maternal Y:
# first, in case of local phasing it is not clear which haplotype is really paternal or maternal
# also, can't interpret het variants seen in vcf-files for Y and MT chrs
# one may want to include these for mapping though
# thus, later the hetsnvs on these chrs should be excluded from read/allele count analyses



$(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_h2.fa: $(OUTPUT_DIR)/maternal.chain
	cat \
		$(OUTPUT_DIR)/1_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/2_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/3_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/4_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/5_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/6_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/7_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/8_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/9_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/10_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/11_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/12_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/13_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/14_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/15_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/16_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/17_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/18_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/19_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/20_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/21_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/22_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/X_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/Y_$(VCF_SAMPLE_ID)_maternal.fa \
		$(OUTPUT_DIR)/MT_$(VCF_SAMPLE_ID)_maternal.fa | \
		sed 's/maternal/h2/g' > $@
 
$(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_h1.fa: $(OUTPUT_DIR)/paternal.chain
	cat \
		$(OUTPUT_DIR)/1_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/2_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/3_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/4_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/5_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/6_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/7_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/8_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/9_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/10_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/11_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/12_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/13_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/14_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/15_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/16_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/17_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/18_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/19_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/20_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/21_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/22_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/X_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/Y_$(VCF_SAMPLE_ID)_paternal.fa \
		$(OUTPUT_DIR)/MT_$(VCF_SAMPLE_ID)_paternal.fa | \
		sed 's/paternal/h1/g' > $@
	

$(OUTPUT_DIR)/bowtie_build_diploid.log: $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_h2.fa $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_h1.fa $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h2.bed $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h1.bed
	$(BOWTIE_build) --offrate 2 \
		$(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_h2.fa,$(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_h1.fa $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_diploid > $@

$(OUTPUT_DIR)/STAR_idx_diploid_Log.out: $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_diploid.$(notdir $(ANNOTATION)) $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_h2.fa $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_h1.fa $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h2.bed $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_h1.bed
	mkdir $(OUTPUT_DIR)/STAR_idx_diploid
	$(STAR) \
		--runThreadN $(N_THREADS) \
		--runMode genomeGenerate \
		--genomeDir $(OUTPUT_DIR)/STAR_idx_diploid \
		--genomeFastaFiles $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_h2.fa $(OUTPUT_DIR)/$(VCF_SAMPLE_ID)_h1.fa \
		--outFileNamePrefix $(OUTPUT_DIR)/STAR_idx_diploid_

