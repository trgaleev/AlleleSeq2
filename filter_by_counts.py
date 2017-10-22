import sys

th_tot  = int(sys.argv[1])
th_each = int(sys.argv[2])
sys.stdout.write(sys.stdin.readline())
for line in sys.stdin:
    (chrm, ref_coord, h1_coord, h2_coord, ref_allele, h1_allele, h2_allele, 
     cA, cC, cG, cT, cN, ref_allele_ratio, sum_ref_n_alt_cnts, p_binom, cnv, mmap_log) = line.split('\t')
    counts = {'A': int(cA), 'C': int(cC), 'G':int(cG), 'T':int(cT), 'N':int(cN)}
    if counts[h1_allele] + counts[h2_allele] >= th_tot and counts[h1_allele] >= th_each and counts[h2_allele] >= th_each:
        sys.stdout.write(line)
