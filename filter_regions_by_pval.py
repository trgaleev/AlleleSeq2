import sys

with open(sys.argv[1], 'r') as infdr:
    for line in infdr:
        if line.startswith('Target'): pth = float(line.split()[-1])

sys.stdout.write(sys.stdin.readline())
for line in sys.stdin:

    region, h1_count, h2_count, h1_allele_ratio, p_binom, snv_count, mmap_log = line.split('\t')
    outcolmns = '\t'.join([region, h1_count, h2_count, h1_allele_ratio, p_binom, snv_count])

    if float(p_binom) <= pth: sys.stdout.write(outcolmns+'\n')
