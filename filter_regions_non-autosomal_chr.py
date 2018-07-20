# the COUNTS_FILE can come with or without non-autosomal chrs (e.g. chrX for female)
# by default, will make sure only autosomal are in; optionally keep ones listed in sys.argv[3:] (e.g. chrX again)

import sys

keep_chrs     = ['chr' + str(c) for c in range(1,23)] + [str(c) for c in range(1,23)]
if len(sys.argv) > 3: 
    keep_chrs = keep_chrs + ['chr' + c for c in sys.argv[3:]] + sys.argv[3:]

region_chr_dict = {}
with open(sys.argv[1],'r') as regionsf:
    for line in regionsf:
        chrm, _, _, region = line.strip().split('\t')
        region_chr_dict [region] = chrm      

with open(sys.argv[2],'w') as rm_regions_f:
    for line in sys.stdin:

        if not line.startswith('#'):
            region, hap1_count, hap2_count, hap1_allele_ratio, p_binom, snv_count = line.split('\t')

            if region_chr_dict[region] in keep_chrs: sys.stdout.write(line)
            else: rm_regions_f.write('\t'.join([region, region_chr_dict[region]] + ['_chr_filter\n']))
        else: sys.stdout.write(line)
