import sys

def check_subset(l1,l2):
    for e in l1:
        if e not in l2: return False
    return True

with open(sys.argv[1], 'a') as rm_hetSNV_f:
    for line in sys.stdin:
        if not line.startswith('#'):
            (chrm,ref_coord,hap1_coord,hap2_coord,ref_allele,hap1_allele,hap2_allele,
                cA,cC,cG,cT,cN,ref_allele_ratio,sum_ref_n_alt_cnts,p_binom,
                warning_hap1,warning_hap2,cnv) = line.split()
            if check_subset(warning_hap1.split(','),['.','zero_cnt']) and check_subset(warning_hap2.split(','),['.','zero_cnt']) and hap1_coord != 'notLifted?' and hap2_coord != 'notLifted?': 
                sys.stdout.write('\t'.join([
                                               chrm,ref_coord,hap1_coord,hap2_coord,ref_allele,hap1_allele,hap2_allele,
                                               cA,cC,cG,cT,cN,ref_allele_ratio,sum_ref_n_alt_cnts,p_binom, cnv
                                           ])+'\n')
            else: 
                rm_hetSNV_f.write('\t'.join([chrm, ref_coord, 
                                    '_'.join([cA,cC,cG,cT,'','um',warning_hap1,warning_hap2])])+'\n')
        else: sys.stdout.write('\t'.join(line.split()[:15] + [line.split()[-1]])+'\n')
