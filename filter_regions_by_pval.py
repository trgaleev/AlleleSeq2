import sys

with open(sys.argv[1], 'r') as infdr:
    for line in infdr:
        if line.startswith('Target'): pth = float(line.split()[-1])

#sys.stdout.write('\t'.join(sys.stdin.readline().split('\t')[0-6]))
for line in sys.stdin:

    region, hap1_count, hap2_count, hap1_allele_ratio, p_binom, snv_count, snv_hap1_hap2_coords, mmap_log = line.split('\t')
    outcolmns = '\t'.join([region, hap1_count, hap2_count, hap1_allele_ratio, p_binom, snv_count])

    if line.startswith('#'): sys.stdout.write(outcolmns+'\n')
    else:
        if float(p_binom) <= pth: sys.stdout.write(outcolmns+'\n')
