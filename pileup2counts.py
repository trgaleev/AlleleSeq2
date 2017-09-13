import binom
import sys
import scipy.stats
import re
import read_pileup
from collections import defaultdict


hetSNV_dict=defaultdict(dict)
hetSNV_list=[]

# read the h1 and h2 coords
with open(sys.argv[2],'r') as in_h1:
    for line in in_h1:
        line = line.strip()
        h_chr, _, h_crd, ref_info = line.split('\t')
        #h_cr, h = h_chr.split('_')
        r_chr, r_crd, r_a, a1, a2 = ref_info.split('_')
        #hetSNV_dict[r_chr+'_'+r_crd][h]={'crd':h_crd}
        hetSNV_dict[r_chr+'_'+r_crd]['h1_pos'] = h_chr+'_'+h_crd
        hetSNV_dict[r_chr+'_'+r_crd]['r_a'] = r_a
        hetSNV_dict[r_chr+'_'+r_crd]['a1'] = a1
        hetSNV_dict[r_chr+'_'+r_crd]['a2'] = a2

with open(sys.argv[3],'r') as in_h2:
    for line in in_h2:
        line=line.strip()
        h_chr, _, h_crd, ref_info = line.split('\t')
        #h_cr, h = h_chr.split('_')
        r_chr, r_crd, r_a, a1, a2 = ref_info.split('_')
        #hetSNV_dict[r_chr+'_'+r_crd][h]={'crd':h_crd}
        hetSNV_dict[r_chr+'_'+r_crd]['h2_pos'] = h_chr+'_'+h_crd
        hetSNV_dict[r_chr+'_'+r_crd]['r_a'] = r_a
        hetSNV_dict[r_chr+'_'+r_crd]['a1'] = a1
        hetSNV_dict[r_chr+'_'+r_crd]['a2'] = a2
        hetSNV_list.append(r_chr+'_'+r_crd)


# read the pileups
pileup_dict = read_pileup.pileup_to_basecnts(sys.argv[4:])

sys.stdout.write('\t'.join([
    '#chr',
    'ref_coord',
    'h1_coord',
    'h2_coord',
    'ref_allele',
    'h1_allele',
    'h2_allele',
    'cA',
    'cC',
    'cG',
    'cT',
    'cN',
    'ref_allele_ratio',
    'sum_ref_n_alt_cnts',
    'p_binom',
    'warn_h1',
    'warn_h2'
])+'\n')

#for k in hetSNV_dict:
# to keep the original order
for k in hetSNV_list:

    basecnts_h1 = pileup_dict.get(hetSNV_dict[k].get('h1_pos',None), {'A':0, 'C':0, 'G':0, 'T':0, 'N':0, 'warning':'zero_cnt'})
    basecnts_h2 = pileup_dict.get(hetSNV_dict[k].get('h2_pos',None), {'A':0, 'C':0, 'G':0, 'T':0, 'N':0, 'warning':'zero_cnt'})
    
    basecnts = {
        'A':basecnts_h1['A'] + basecnts_h2['A'],
        'C':basecnts_h1['C'] + basecnts_h2['C'],
        'G':basecnts_h1['G'] + basecnts_h2['G'],
        'T':basecnts_h1['T'] + basecnts_h2['T'],
        'N':basecnts_h1['N'] + basecnts_h2['N']
    }

    tot_cnt = basecnts[hetSNV_dict[k]['a1']]+basecnts[hetSNV_dict[k]['a2']]

    if tot_cnt >= int(sys.argv[1]):

        ref_cnt = basecnts[hetSNV_dict[k]['r_a']]
	if ref_cnt > tot_cnt:
		sys.stderr.write(sys.argv[0] +' WARNING:\t' + k + ', ref allele cnt (' + str(ref_cnt) + ') is larger than ref + alt allele counts (' + str(tot_cnt) + ')\n')
		sys.stderr.write(str(basecnts) + '\n')
		sys.stderr.write('miscalled multi-allelic variant? Skipping this hetSNV\n')
		continue
	#print basecnts
	#print k +'\t' + str(ref_cnt) + '\t' + str(tot_cnt) + '\n'
        pbinom=binom.binomtest(ref_cnt, tot_cnt, 0.5)
        sys.stdout.write('\t'.join([    
            k.split('_')[0],
            k.split('_')[1],
            hetSNV_dict[k].get('h1_pos','notLifted?'),
            hetSNV_dict[k].get('h2_pos','notLifted?'),
            hetSNV_dict[k]['r_a'],
            hetSNV_dict[k]['a1'],
            hetSNV_dict[k]['a2'],
            str(basecnts['A']),        
            str(basecnts['C']),        
            str(basecnts['G']),        
            str(basecnts['T']),
            str(basecnts['N']),
            str(float(ref_cnt)/float(tot_cnt)),
            str(tot_cnt),
            str(pbinom),
            basecnts_h1['warning'],
            basecnts_h2['warning']
        ])+'\n')

