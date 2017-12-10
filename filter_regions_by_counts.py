import sys

th_tot  = int(sys.argv[1])
th_each = int(sys.argv[2])
sys.stdout.write(sys.stdin.readline())
for line in sys.stdin:

    region,h1_count,h2_count,h1_allele_ratio,p_binom,snv_count,mmap_log = line.split('\t')
    if int(h1_count) + int(h2_count) >= th_tot and int(h1_count) >= th_each and int(h2_count) >= th_each:
        sys.stdout.write(line)
