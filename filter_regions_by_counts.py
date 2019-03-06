import sys

th_tot  = int(sys.argv[1])
th_each = int(sys.argv[2])
sys.stdout.write(sys.stdin.readline())
for line in sys.stdin:

    region,hap1_count,hap2_count,hap1_allele_ratio,p_binom,snv_count,snv_hap1_hap2_coords,mmap_log = line.split('\t')
    if int(hap1_count) + int(hap2_count) >= th_tot and int(hap1_count) >= th_each and int(hap2_count) >= th_each:
        sys.stdout.write(line)
