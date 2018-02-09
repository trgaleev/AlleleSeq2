import sys

e_d = {}
for line in sys.stdin:
    chrm, _, crd, e = line.split()
    e_d[chrm+'_'+crd] = e


h1f = open(sys.argv[2],'w')
h2f = open(sys.argv[3],'w')

with open(sys.argv[1],'r') as cntf:
    cntf.readline()
    for line in cntf:
        (chrm,ref_coord,h1_coord,h2_coord,ref_allele,h1_allele,h2_allele,
        cA,cC,cG,cT,cN,ref_allele_ratio,sum_ref_n_alt_cnts,p_binom,cnv,mmap_log) = line.split()

        if chrm+'_'+ref_coord in e_d:

            h1_coord = h1_coord.split('_')[-1]
            h2_coord = h2_coord.split('_')[-1]

            h1f.write('\t'.join([chrm+'_h1', str(int(h1_coord)-1), h1_coord, e_d[chrm+'_'+ref_coord]])+'\n')
            h2f.write('\t'.join([chrm+'_h2', str(int(h2_coord)-1), h2_coord, e_d[chrm+'_'+ref_coord]])+'\n')

h1f.close()
h2f.close()

