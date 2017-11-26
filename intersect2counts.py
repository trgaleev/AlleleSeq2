import binom
import sys
import scipy.stats


e_list = []
e_dict = {}
for line in sys.stdin:
    #chr1_h1    1055394    1055395    D00777:83:C8R9PACXX:5:2107:17632:19962/1    255    -    1055295    1055396    0,0,0    1    101,    0,    1    chr1_h1    1055394    1055395    ENSG00000188157.13
    # chr1_h1	909877	909878	chr1_909764_910952	5	chr1_h1	909833	909909	GLP3B7-C38:1381:HY5WKBCXX:2:1108:3941:50109/1	255	+

    rname = line.split('\t')[8].split('/')[0]
    ename = line.split('\t')[3]
    hap = line.split('\t')[0].split('_')[1]
    snv = line.split('\t')[0]+'_'+line.split()[2]


    if ename not in e_dict:
        e_list.append(ename)
        e_dict[ename]={'h1':set(), 'h2':set(), 'snv_count': 0, 'snvs':':'}

    e_dict[ename][hap].add(rname)
    if snv not in e_dict[ename]['snvs']: e_dict[ename]['snvs'] += snv + ':'


with open('hets_intervals_h1.bed', 'r') as in1:
    for line in in1:
        if line.split()[-1] in e_dict: e_dict[line.split()[-1]]['snv_count'] += 1


if sys.argv[1] == 'uniq':
    sys.stdout.write('\t'.join(['#chr',
                                'start_coord',
                                'end_coord',
                                'h1_count',
                                'h2_count',
                                'h1_allele_ratio',
                                'p_binom',
                                'snv_count',
                                'snv_h1_h2_coords']) + '\n')

    for e in e_list:
        h1_count = len(e_dict[e]['h1'])
        h2_count = len(e_dict[e]['h2'])

        pbinom=binom.binomtest(h1_count, h1_count+h2_count, 0.5)

        h1_allele_ratio = float(h1_count)/(float(h1_count)+float(h2_count))

        sys.stdout.write('\t'.join([
                    '\t'.join(e.split('_')),
                    str(h1_count),
                    str(h2_count),
                    str(h1_allele_ratio),
                    str(pbinom),
                    str(e_dict[e]['snv_count']),
                    e_dict[e]['snvs']
                    ])+'\n')

    

elif sys.argv[1] == 'mmap':
    sys.stdout.write('\t'.join(['#chr',
                                'start_coord',
                                'end_coord',
                                'h1_count',
                                'h2_count']) + '\n')
    for e in e_list:
        h1_count = len(e_dict[e]['h1'])
        h2_count = len(e_dict[e]['h2'])


        sys.stdout.write('\t'.join([
                    '\t'.join(e.split('_')),
                    str(h1_count),
                    str(h2_count)
                    ])+'\n')



