import sys
import numpy

snv_list = []
snv_dict = {}
for line in sys.stdin:
    chr, pos0, pos1, inf, cnt = line.split()
    snv_list.append('\t'.join([inf]))
    snv_dict['\t'.join([inf])] = int(cnt)

rd_median = numpy.median(snv_dict.values())

sys.stdout.write('#chr\tpos\trd\n')
for snv in snv_list:
    sys.stdout.write('\t'.join(snv.split('_')[:2]+[str(snv_dict[snv]/rd_median)])+'\n')

sys.stderr.write('rd_median:\t'+str(rd_median)+'\n')
