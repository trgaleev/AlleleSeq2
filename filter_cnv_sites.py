import sys

cnv_dict={}
with open(sys.argv[2], 'r') as cnvfile:
    for line in cnvfile:
        if not line.startswith('chrm'):
            chrm,snppos,rd = line.strip().split('\t')
            cnv_dict[chrm+'_'+snppos] = rd


with open(sys.argv[1], 'a') as rm_hetSNV_f:
    for line in sys.stdin:
        if not line.startswith('#'):
            (chrm,ref_coord,hap1_coord,hap2_coord,ref_allele,hap1_allele,hap2_allele,
                cA,cC,cG,cT,cN,ref_allele_ratio,sum_ref_n_alt_cnts,p_binom,
                warning_hap1,warning_hap2) = line.split()
            rd = float(cnv_dict[chrm+"_"+ref_coord])
            if 0.5 <= rd <= 1.5:
                sys.stdout.write(line.strip()+'\t'+str(rd)+'\n')
            else: 
                rm_hetSNV_f.write('\t'.join([chrm, ref_coord, 
                                     '_'.join([cA,cC,cG,cT,'','rd:'+str(rd)])])+'\n')
        else: sys.stdout.write(line.strip()+'\tcnv\n')
