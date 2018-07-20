import sys

th_tot  = int(sys.argv[1])
th_each = int(sys.argv[2])
sys.stdout.write(sys.stdin.readline())
for line in sys.stdin:
    (chrm, ref_coord, hap1_coord, hap2_coord, ref_allele, hap1_allele, hap2_allele, 
     cA, cC, cG, cT, cN, ref_allele_ratio, sum_ref_n_alt_cnts, p_binom, cnv, mmap_log) = line.split('\t')
    counts = {'A': int(cA), 'C': int(cC), 'G':int(cG), 'T':int(cT), 'N':int(cN)}
    if counts[hap1_allele] + counts[hap2_allele] >= th_tot and counts[hap1_allele] >= th_each and counts[hap2_allele] >= th_each:
        sys.stdout.write(line)
