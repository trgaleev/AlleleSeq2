import sys

e_d = {}
for line in sys.stdin:
    chrm, _, crd, e = line.split()
    e_d[chrm+'_'+crd] = e


hap1f = open(sys.argv[2],'w')
hap2f = open(sys.argv[3],'w')

with open(sys.argv[1],'r') as cntf:
    cntf.readline()
    for line in cntf:
        (chrm,ref_coord,hap1_coord,hap2_coord,ref_allele,hap1_allele,hap2_allele,
        cA,cC,cG,cT,cN,ref_allele_ratio,sum_ref_n_alt_cnts,p_binom,cnv,mmap_log) = line.split()

        if chrm+'_'+ref_coord in e_d:

            hap1_coord = hap1_coord.split('_')[-1]
            hap2_coord = hap2_coord.split('_')[-1]

            hap1f.write('\t'.join([chrm+'_hap1', str(int(hap1_coord)-1), hap1_coord, e_d[chrm+'_'+ref_coord]])+'\n')
            hap2f.write('\t'.join([chrm+'_hap2', str(int(hap2_coord)-1), hap2_coord, e_d[chrm+'_'+ref_coord]])+'\n')

hap1f.close()
hap2f.close()

