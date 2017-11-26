import sys

with open(sys.argv[1], 'r') as infdr:
    for line in infdr:
        if line.startswith('Target'): pth = float(line.split()[-1])

for line in sys.stdin:
    (chr,ref_coord,h1_coord,h2_coord,ref_allele,h1_allele,h2_allele,cA,cC,cG,cT,cN,
    ref_allele_ratio,sum_ref_n_alt_cnts,p_binom,cnv,mmap_log) = line.split('\t')
    outcolmns = '\t'.join([
                                       chr,
                                       ref_coord,
                                       ref_allele,
                                       h1_allele,
                                       h2_allele,
                                       cA, cC, cG, cT, cN,
                                       ref_allele_ratio,
                                       p_binom
                                    ])

    if chr == '#chr': sys.stdout.write(outcolmns+'\n')
    elif float(p_binom) <= pth: sys.stdout.write(outcolmns+'\n')
