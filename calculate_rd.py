import sys
import numpy
from collections import defaultdict

snv_dict = defaultdict(int)

for line in sys.stdin:
    inf = line.split()[-1]
    snv_dict[inf] += 1

rd_median = numpy.median(snv_dict.values())

sys.stdout.write('#chr\tpos\trd\n')
for snv in snv_dict:
    sys.stdout.write('\t'.join(snv.split('_')[:2]+[str(snv_dict[snv]/rd_median)])+'\n')

sys.stderr.write('median rd:\t'+str(rd_median)+'\n')
