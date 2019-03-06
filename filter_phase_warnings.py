import sys

def check_subset(l1,l2):
    for e in l1:
        if e not in l2: return False
    return True

with open(sys.argv[1], 'w') as rm_hetSNV_f:
#rm_hetSNV_f = open(sys.argv[1], 'a') 
    rm_hetSNV_f.write('#chr\tref_coord\tcA_cC_cG_cT__warning_hap1__warning_hap2\n')
    for line in sys.stdin:
        if not line.startswith('#'):
            (chrm,ref_coord,hap1_coord,hap2_coord,ref_allele,hap1_allele,hap2_allele,
                cA,cC,cG,cT,cN,ref_allele_ratio,sum_ref_n_alt_cnts,p_binom,
                warning_hap1,warning_hap2,cnv) = line.split()
	    # not using 'zero_cnt' as a warning anymore
            # since there is a custom filtering for a min of reads supporting each allele further in the pipeline
            #if check_subset(warning_hap1.split(','),['.','zero_cnt']) and check_subset(warning_hap2.split(','),['.','zero_cnt']) and hap1_coord != 'not_lifted?' and hap2_coord != 'not_lifted?': 
            if check_subset(warning_hap1.split(','),['.']) and check_subset(warning_hap2.split(','),['.']) and hap1_coord != 'not_lifted?' and hap2_coord != 'not_lifted?': 
                sys.stdout.write('\t'.join([
                                               chrm,ref_coord,hap1_coord,hap2_coord,ref_allele,hap1_allele,hap2_allele,
                                               cA,cC,cG,cT,cN,ref_allele_ratio,sum_ref_n_alt_cnts,p_binom, cnv
                                           ])+'\n')
            else: 
                rm_hetSNV_f.write('\t'.join([chrm, ref_coord, 
                                    '_'.join([cA,cC,cG,cT,'',warning_hap1,'',warning_hap2])])+'\n')
        else: sys.stdout.write('\t'.join(line.split()[:15] + [line.split()[-1]])+'\n')

#rm_hetSNV_f.close()
