import binom
import sys
import scipy.stats


e_list = []
e_dict = {}
for line in sys.stdin:
	#chr1_h1	1055394	1055395	D00777:83:C8R9PACXX:5:2107:17632:19962/1	255	-	1055295	1055396	0,0,0	1	101,	0,	1	chr1_h1	1055394	1055395	ENSG00000188157.13
	rname = line.split()[7].split('/')[0]
	ename = line.split()[3]
	hap = line.split()[0].split('_')[1]
	snv = line.split()[0]+'_'+line.split()[1]
	if ename not in e_dict:
		e_list.append(ename)
		e_dict[ename]={'h1':{}, 'h2':{}, 'snv_count': 0, 'snvs':':'}

	e_dict[ename][hap][rname] = None
	e_dict[ename]['snvs'] = e_dict[ename]['snvs'] + snv + ':'


with open('hets_intervals_h1.bed', 'r') as in1:
	for line in in1:
		e_dict[line.split()[-1]]['snv_count'] += 1

sys.stdout.write('region\th1_count\th2_count\th1_allele_ratio\tp_binom\tsnv_count\tsnv_coords\n')
for e in e_list:
	h1_count = len(e_dict[e]['h1'])
	h2_count = len(e_dict[e]['h2'])

	pbinom=binom.binomtest(h1_count, h1_count+h2_count, 0.5)

	h1_allele_ratio = float(h1_count)/(float(h1_count)+float(h2_count))

	sys.stdout.write('\t'.join([
					e,
					str(h1_count),
					str(h2_count),
					str(h1_allele_ratio),
					str(pbinom),
					str(e_dict[e]['snv_count']),
					e_dict[e]['snvs']

			
					])+'\n')

	

