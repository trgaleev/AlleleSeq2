The ENTEx project generates bam files of RNA-seq and ChIP-seq where the reads are aligned to both haplotypes of the personal genomes. The script calc_read_coverage.sh calculates the read coverage on each haplotype and converts the coordinates in the personal genomes to the reference genome (GRCh38). It also generate average read coverage of bins of 2.5 Mb, which is used to make the chromosome painting (see git repo gersteinlab/ChromosomePaintingTool).

Users shoud have samtools (v1.9 was tested), bedtools (v 2.27.1 was tested), liftOver (available in UCSC genome browser utitlites), and bedGraphToBigWig (available in UCSC genome browser utitlites) installed and added to the system path. 

The "data" directory contains chain files used to convert the coordinates in the personal genomes to the reference genome for the 4 ENTEx individuals, and other necessary files.

To use, download all content and give permission to run calc_read_coverage.sh. On Linux, permission is granted by the following command:

```bash
chmod u+x calc_read_coverage.sh
```
**To run the script, provide 5 arguemnts:**
```bash
./calc_read_coverage.sh path_to_bam_files path_to_output prefix_of_output_files path_to_chain_files N_cpus
```
If there are multiple bam files (e.g. replicates of the same assay), the script will calculate a total read coverage from all the bam files. Merging bam files can benefit from using multiple CPUs, but is not required.


**Output:**

prefix.sorted.merged_hap_uniq.bedgraph --- this is the read coverage on both haplotypes, using coordinates in the personal genome

prefix.sorted.hapXonref.bedgraph --- read coverage of haplotype X, using coordinates in the reference genome

prefix.sorted.hapXonref.bigwig --- same as above, but in bigwig format

prefix.binned.hapX.bedgraph --- average read coverage of bins of 2.5 Mb of haplotype X, using coordinates in the reference genome

prefix.Haplotype1.txt --- file used by the chromosome painting tool