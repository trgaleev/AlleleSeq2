import sys

if 'interval_FDR.binom.txt' in sys.argv[1]:
    with open(sys.argv[1], 'r') as infdr:
        for line in infdr:
            if line.startswith('Target'): pth = float(line.split()[-1])
else: # must be a number
    pth = float(sys.argv[1])


sys.stdin.readline()
for line in sys.stdin:

    region,h1_count,h2_count,h1_allele_ratio,p_binom,snv_count,snv_coords = line.split('\t')
    outcolmns = '\t'.join([region,h1_count,h2_count,h1_allele_ratio,p_binom,snv_count,snv_coords])

    if float(p_binom) <= pth: sys.stdout.write(outcolmns+'\n')
