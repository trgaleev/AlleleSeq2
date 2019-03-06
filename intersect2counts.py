import binom
import sys
import scipy.stats


e_list = []
e_dict = {}
for line in sys.stdin:


    # seems that, for single -b bam file, the formatting is:
    # chr1_hap1	187152	187153	ENSG00000279457.3	chr1_hap1	186410	187170	D00777:64:C84FFACXX:4:2109:3419:88661/1	255	+
    # for two bam files (the new colmn shows which file?)
    # chr1_hap1	912441	912442	chr1_912852_913352	2	chr1_hap1	912401	912477	GLP0C0-857:1553:HHVYFBCXY:2:1115:1549:46799	255	+

    if    len(line.split('\t')) == 10: rname = line.split('\t')[7].split('/')[0]   # the last split is to not double-count paired ends of the same read
    elif  len(line.split('\t')) == 11: rname = line.split('\t')[8].split('/')[0]  
    else: sys.exit(sys.argv[0] + ':\t confused with intersecBed output:\n' + line)  

    ename = line.split('\t')[3]
    hap = line.split('\t')[0].split('_')[1]
    snv = line.split('\t')[0] + '_' + line.split()[2]

    if ename not in e_dict:
        e_list.append(ename)
        e_dict[ename]={'hap1': set(), 'hap2': set(), 'snv_count': 0, 'snvs': ':'}

    # todo: may not include either hap's coord, if no reads mapped to it
    # so, for some regions in the 'snv_hap1_hap2_coords', some hetSNV may have coord from only one hap 
    # remove this column here in all future steps if it isn't necessary, or get this info from the bed file

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



