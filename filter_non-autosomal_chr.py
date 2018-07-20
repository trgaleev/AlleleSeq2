import sys

keep_chrs     = ['chr' + str(c) for c in range(1,23)] + [str(c) for c in range(1,23)]
if len(sys.argv) > 2: 
    keep_chrs = keep_chrs + ['chr' + c for c in sys.argv[2:]] + sys.argv[2:]


with open(sys.argv[1],'a') as rm_hetSNV_f:
    for line in sys.stdin:

        if not line.startswith('#'):
            (chrm,ref_coord,hap1_coord,hap2_coord,ref_allele,hap1_allele,hap2_allele,
                     cA,cC,cG,cT,cN,ref_allele_ratio,sum_ref_n_alt_cnts,p_binom,
                     warning_hap1,warning_hap2,cnv) = line.split('\t')

            if chrm in keep_chrs: sys.stdout.write(line)
            else: rm_hetSNV_f.write('\t'.join([chrm,ref_coord, 
                   '_'.join([cA,cC,cG,cT,'','chr_filter'])])+'\n')
        else: sys.stdout.write(line)
