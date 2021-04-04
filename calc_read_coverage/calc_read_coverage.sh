#!/bin/bash

#make temporary output directory
mkdir "$2"/tmp

#################################################################################
### preproces bam files
#################################################################################
echo processing bam...

# exclude multi-mapping reads
for i in $(ls "$1"/*.bam)
do
	filename=$(echo $i | sed 's/\.bam//g')
	samtools view -@ $5 -h -q 255 -b $i > "$filename".hap_uniq.bam	
done

# if there is only one bam, index it
N=$(ls "$1"/*.hap_uniq.bam | wc -l)
if ((N==1)); then
    file=$(ls "$1"/*hap_uniq.bam)
    samtools index -@ $5 $file  
else
#merge bam files and index
    samtools merge -@ $5 "$1"/merged_hap_uniq.bam "$1"/*hap_uniq.bam
    samtools index -@ $5 "$1"/merged_hap_uniq.bam
    file=$(echo "$1"/merged_hap_uniq.bam)
fi

#################################################################################
### calculate coverage from bam
#################################################################################
echo calculating coverage...
bedtools genomecov -ibam $file -bg -split > "$2"/tmp/merged_hap_uniq.bedgraph

#################################################################################
### generate bedgraph and bigwig that mapped to the reference genome
#################################################################################

#renaming chrMT to chrM
sed -i 's/chrMT/chrM/g' "$2"/tmp/merged_hap_uniq.bedgraph

#liftOver to reference. hapXtoref.chain is the same as ${individual}.hap1.swaped.chain
echo lifting over to reference...
liftOver "$2"/tmp/merged_hap_uniq.bedgraph "$4"/hap1toref.chain "$2"/"$3".hap1onref.bedgraph "$2"/tmp/"$3".unmapped.hap1
liftOver "$2"/tmp/merged_hap_uniq.bedgraph "$4"/hap2toref.chain "$2"/"$3".hap2onref.bedgraph "$2"/tmp/"$3".unmapped.hap2

#sort bedgraphs
echo sorting...
bedtools sort -i "$2"/tmp/merged_hap_uniq.bedgraph > "$2"/"$3".sorted.merged_hap_uniq.bedgraph
bedtools sort -i "$2"/"$3".hap1onref.bedgraph > "$2"/"$3".sorted.hap1onref.bedgraph
bedtools sort -i "$2"/"$3".hap2onref.bedgraph > "$2"/"$3".sorted.hap2onref.bedgraph

#convert bedgraph to bigwig
echo making bigwig...
bedGraphToBigWig "$2"/"$3".sorted.hap1onref.bedgraph "$4"/GRCh38_chrom_sizes.genome "$2"/"$3".hap1onref.bigwig
bedGraphToBigWig "$2"/"$3".sorted.hap2onref.bedgraph "$4"/GRCh38_chrom_sizes.genome "$2"/"$3".hap2onref.bigwig

#clean up
rm -r "$2"/tmp
mv "$2"/"$3".sorted.hap1onref.bedgraph "$2"/"$3".hap1onref.bedgraph
mv "$2"/"$3".sorted.hap2onref.bedgraph "$2"/"$3".hap2onref.bedgraph

#################################################################################
### the code below generate file for chromosome painting 
#################################################################################

# calcuate average coverage in 2.5Mb bins
echo binnning...
bedtools map -a "$4"/25Mb_hg38_bins.bed -b "$2"/"$3".hap1onref.bedgraph -c 4 -o mean \
> "$2"/"$3".mapped.hap1.bedgraph
bedtools map -a "$4"/25Mb_hg38_bins.bed -b "$2"/"$3".hap2onref.bedgraph -c 4 -o mean \
> "$2"/"$3".mapped.hap2.bedgraph

#generate files for chromosome painting 
python annotation_file_maker.1.py "$2"/"$3".mapped.hap1.bedgraph "$2"/"$3".mapped.hap2.bedgraph "$2"/"$3".Haplotype1.txt "$2"/"$3".Haplotype2.txt
python annotation_file_maker.2.py "$2"/"$3".mapped.hap1.bedgraph "$2"/"$3".mapped.hap2.bedgraph "$2"/"$3".Haplotype1.txt "$2"/"$3".Haplotype2.txt

#make a copy of the bedgraph for the paingting tool, which requires a different file name
cp "$2"/"$3".hap1onref.bedgraph "$2"/"$3".Haplotype1.bedgraph
cp "$2"/"$3".hap2onref.bedgraph "$2"/"$3".Haplotype2.bedgraph

#remove tmp files
rm "$2"/"$3".mapped.hap*.bedgraph "$2"/"$3".mapped.hap*.bedgraph_list

echo done
