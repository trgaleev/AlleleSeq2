import sys
with open(sys.argv[1],'a') as rm_hetSNV_f:
    for line in sys.stdin:
        if not line.startswith('#'):
            if not any(x in line.split()[0] for x in ['X','Y','M']): 
                sys.stdout.write(line)
            else: 
                (chr,ref_coord,h1_coord,h2_coord,ref_allele,h1_allele,h2_allele,
                     cA,cC,cG,cT,cN,ref_allele_ratio,sum_ref_n_alt_cnts,p_binom,
                     warning_h1,warning_h2,cnv) = line.split('\t')
                rm_hetSNV_f.write('\t'.join([chr,ref_coord, 
                                  '_'.join([cA,cC,cG,cT,'','chr_filter'])])+'\n')
        else: sys.stdout.write(line)
