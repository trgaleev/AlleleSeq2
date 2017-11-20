import sys

def check_subset(l1,l2):
    for e in l1:
        if e not in l2: return False
    return True

with open(sys.argv[1], 'a') as rm_hetSNV_f:
    for line in sys.stdin:
        if not line.startswith('#'):
            (chrm,ref_coord,h1_coord,h2_coord,ref_allele,h1_allele,h2_allele,
                cA,cC,cG,cT,cN,ref_allele_ratio,sum_ref_n_alt_cnts,p_binom,
                warning_h1,warning_h2,cnv) = line.split()
            if check_subset(warning_h1.split(','),['.','zero_cnt']) and check_subset(warning_h2.split(','),['.','zero_cnt']) and h1_coord != 'notLifted?' and h2_coord != 'notLifted?': 
                sys.stdout.write('\t'.join([
                                               chrm,ref_coord,h1_coord,h2_coord,ref_allele,h1_allele,h2_allele,
                                               cA,cC,cG,cT,cN,ref_allele_ratio,sum_ref_n_alt_cnts,p_binom, cnv
                                           ])+'\n')
            else: 
                rm_hetSNV_f.write('\t'.join([chrm, ref_coord, 
                                    '_'.join([cA,cC,cG,cT,'','um',warning_h1,warning_h2])])+'\n')
        else: sys.stdout.write('\t'.join(line.split()[:15] + [line.split()[-1]])+'\n')
