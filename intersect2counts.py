import binom
import sys
import scipy.stats


e_list = []
e_dict = {}
for line in sys.stdin:
    #chr1_hap1    1055394    1055395    D00777:83:C8R9PACXX:5:2107:17632:19962/1    255    -    1055295    1055396    0,0,0    1    101,    0,    1    chr1_hap1    1055394    1055395    ENSG00000188157.13
    # chr1_hap1	909877	909878	chr1_909764_910952	5	chr1_hap1	909833	909909	GLP3B7-C38:1381:HY5WKBCXX:2:1108:3941:50109/1	255	+


    # chr1_hap1	187152	187153	ENSG00000279457.3	chr1_hap1	186410	187170	D00777:64:C84FFACXX:4:2109:3419:88661/1	255	+

    rname = line.split('\t')[7].split('/')[0]   # the last split is to not double-count paired ends of the same read
    ename = line.split('\t')[3]
    hap = line.split('\t')[0].split('_')[1]
    snv = line.split('\t')[0] + '_' + line.split()[2]

    if ename not in e_dict:
        e_list.append(ename)
        e_dict[ename]={'hap1': set(), 'hap2': set(), 'snv_count': 0, 'snvs': ':'}

    # may not include either hap's coord, if no reads mapped to it
    # todo: remove if this isn't necessary, or get this info from the bed file

    e_dict[ename][hap].add(rname)
    if snv not in e_dict[ename]['snvs']: e_dict[ename]['snvs'] += snv + ':'


with open(sys.argv[2], 'r') as in1:
    for line in in1:
        if line.split()[-1] in e_dict: e_dict[line.split()[-1]]['snv_count'] += 1


if sys.argv[1] == 'uniq':
    sys.stdout.write('\t'.join(['#region',
                                'hap1_count',
                                'hap2_count',
                                'hap1_allele_ratio',
                                'p_binom',
                                'snv_count',
                                'snv_hap1_hap2_coords']) + '\n')

    for e in e_list:
        hap1_count = len(e_dict[e]['hap1'])
        hap2_count = len(e_dict[e]['hap2'])

        pbinom=binom.binomtest(hap1_count, hap1_count+hap2_count, 0.5)

        hap1_allele_ratio = float(hap1_count)/(float(hap1_count)+float(hap2_count))

        sys.stdout.write('\t'.join([
                    e,
                    str(hap1_count),
                    str(hap2_count),
                    str(hap1_allele_ratio),
                    str(pbinom),
                    str(e_dict[e]['snv_count']),
                    e_dict[e]['snvs']
                    ])+'\n')
    

elif sys.argv[1] == 'mmap':
    sys.stdout.write('\t'.join(['#region',
                                'hap1_count',
                                'hap2_count']) + '\n')
    for e in e_list:
        hap1_count = len(e_dict[e]['hap1'])
        hap2_count = len(e_dict[e]['hap2'])


        sys.stdout.write('\t'.join([
                    e,
                    str(hap1_count),
                    str(hap2_count)
                    ])+'\n')



