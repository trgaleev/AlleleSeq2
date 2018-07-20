import binom
import sys
import scipy.stats
import re
import read_pileup
from collections import defaultdict


hetSNV_dict = defaultdict(dict)
bases = set(['A', 'C', 'G', 'T', 'N'])
discarded = set([])

rmvdhets_file = open(sys.argv[5],'w')

duplicate_posns = {}
# read the hap1 and hap2 coords
with open(sys.argv[3],'r') as in_hap1:
    for line in in_hap1:
        line = line.strip()
        h_chr, _, h_crd, ref_info = line.split('\t')
        #h_cr, h = h_chr.split('_')
        r_chr, r_crd, r_a, a1, a2 = ref_info.split('_')
        if a1 not in bases or a2 not in bases: 
            rmvdhets_file.write('\t'.join([r_chr, r_crd, ref_info]) + '\n')
            discarded.add(r_chr + '_' + r_crd)
            continue
        if 'hap1_pos' in hetSNV_dict[r_chr+'_'+r_crd] or r_chr+'_'+r_crd in duplicate_posns:
            duplicate_posns[r_chr+'_'+r_crd] = None
            continue
        hetSNV_dict[r_chr+'_'+r_crd]['hap1_pos'] = h_chr + '_' + h_crd
        hetSNV_dict[r_chr+'_'+r_crd]['r_a'] = r_a
        hetSNV_dict[r_chr+'_'+r_crd]['hap1_a'] = a1
        hetSNV_dict[r_chr+'_'+r_crd]['hap2_a'] = a2
with open(sys.argv[4],'r') as in_hap2:
    for line in in_hap2:
        line = line.strip()
        h_chr, _, h_crd, ref_info = line.split('\t')
        #h_cr, h = h_chr.split('_')
        r_chr, r_crd, r_a, a1, a2 = ref_info.split('_')
        if a1 not in bases or a2 not in bases: 
            if r_chr + '_' + r_crd not in discarded: rmvdhets_file.write('\t'.join([r_chr, r_crd, ref_info]) + '\n')
            continue
        if 'hap2_pos' in hetSNV_dict[r_chr+'_'+r_crd] or r_chr+'_'+r_crd in duplicate_posns:
            duplicate_posns[r_chr+'_'+r_crd] = None
            continue
        hetSNV_dict[r_chr+'_'+r_crd]['hap2_pos'] = h_chr + '_' + h_crd
        hetSNV_dict[r_chr+'_'+r_crd]['r_a'] = r_a
        hetSNV_dict[r_chr+'_'+r_crd]['hap1_a'] = a1  # mostly redundant, but needed for cases, where
        hetSNV_dict[r_chr+'_'+r_crd]['hap2_a'] = a2  # an snv couldn't be lifted over for one of the haplotypes

for i in duplicate_posns:
    del hetSNV_dict[i]
    sys.stderr.write('WARN: same coords of two SNVs: ' + ' '.join(i.split('_')) + ' ? Removing.\n')
    rmvdhets_file.write('\t'.join(i.split('_') + ['duplicate_pos']) + '\n')
    # this doesn't catch homozygous SNVs though, when duplicated with a het - vcf2diploid may actually incorporate either? version
    # thus, the input vcf should not have duplicated snvs
    # todo: add a check when generating ref.bed

# read the pileups
pileup_dict = read_pileup.pileup_to_basecnts(sys.argv[6:])

sys.stdout.write('\t'.join([
    '#chr',
    'ref_coord',
    'hap1_coord',
    'hap2_coord',
    'ref_allele',
    'hap1_allele',
    'hap2_allele',
    'cA',
    'cC',
    'cG',
    'cT',
    'cN',
    'ref_allele_ratio',
    'sum_ref_n_alt_cnts',
    'p_binom',
    'warn_hap1',
    'warn_hap2'
])+'\n')


with open(sys.argv[2], 'r') as in_ref:
    for line in in_ref:

        k = line.split('\t')[0] + '_' + line.split('\t')[2]

        # when not lifted neither to hap1 nor hap2
        # or a non A, C, G, T, N nt in vcf
        if k not in hetSNV_dict: continue

        basecnts_hap1 = pileup_dict.get(hetSNV_dict[k].get('hap1_pos',None), {'A':0, 'C':0, 'G':0, 'T':0, 'N':0, 'warning':'zero_cnt', 'hap1_a':'.'})
        basecnts_hap2 = pileup_dict.get(hetSNV_dict[k].get('hap2_pos',None), {'A':0, 'C':0, 'G':0, 'T':0, 'N':0, 'warning':'zero_cnt', 'hap2_a':'.'})

        basecnts = {
            'A':basecnts_hap1['A'] + basecnts_hap2['A'],
            'C':basecnts_hap1['C'] + basecnts_hap2['C'],
            'G':basecnts_hap1['G'] + basecnts_hap2['G'],
            'T':basecnts_hap1['T'] + basecnts_hap2['T'],
            'N':basecnts_hap1['N'] + basecnts_hap2['N']
        }

        #tot_cnt = basecnts.get(basecnts_hap1['hap1_a'], 0) + basecnts.get(basecnts_hap2['hap2_a'], 0)
        tot_cnt =  basecnts[hetSNV_dict[k]['hap1_a']] +  basecnts[hetSNV_dict[k]['hap2_a']]

        if tot_cnt >= max(int(sys.argv[1]), 1):
        # need at least one read to get at least one hap nt from mpileups
        # using that to figure out what phasing was assigned by vcf2diploid
        # if wasn't phased in the .vcf:

            if   (basecnts_hap1['hap1_a'], basecnts_hap2['hap2_a']) in [(hetSNV_dict[k]['hap1_a'], hetSNV_dict[k]['hap2_a']), (hetSNV_dict[k]['hap1_a'], '.'), ('.', hetSNV_dict[k]['hap2_a'])]:
                basecnts_hap1['hap1_a'], basecnts_hap2['hap2_a'] = hetSNV_dict[k]['hap1_a'], hetSNV_dict[k]['hap2_a']
                
            elif (basecnts_hap2['hap2_a'], basecnts_hap1['hap1_a']) in [(hetSNV_dict[k]['hap1_a'], hetSNV_dict[k]['hap2_a']), (hetSNV_dict[k]['hap1_a'], '.'), ('.', hetSNV_dict[k]['hap2_a'])]:
                basecnts_hap2['hap2_a'], basecnts_hap1['hap1_a'] = hetSNV_dict[k]['hap1_a'], hetSNV_dict[k]['hap2_a']

            else:
		#sys.exit(sys.argv[0] + '\n confused with phasing in vcf vs mpileup / hap sequences\n' + line)
		sys.stderr.write(sys.argv[0] + '   WARN: confused with phasing in .vcf vs .mpileup, skipping: ' + k + '\n')
                rmvdhets_file.write('\t'.join(k.split('_') + ['vcf_vs_mpileup_alleles']) + '\n')
                continue


            ref_cnt = basecnts[hetSNV_dict[k]['r_a']]
            if ref_cnt > tot_cnt:  # miscalled multi-allelic variant? 
                rmvdhets_file.write(k.split('_')[0] + '\t' + k.split('_')[1] + '\tref_allele:' + str(ref_cnt) + '_hap1:' + str(basecnts[hetSNV_dict[k]['hap1_a']]) + '_hap2:' + str(basecnts[hetSNV_dict[k]['hap2_a']]) + '\n')
                continue

            pbinom = binom.binomtest(ref_cnt, tot_cnt, 0.5)
            sys.stdout.write('\t'.join([    
                k.split('_')[0],
                k.split('_')[1],
                hetSNV_dict[k].get('hap1_pos','notLifted?'),
                hetSNV_dict[k].get('hap2_pos','notLifted?'),
                hetSNV_dict[k]['r_a'],
                basecnts_hap1['hap1_a'],
                basecnts_hap2['hap2_a'],
                str(basecnts['A']),        
                str(basecnts['C']),        
                str(basecnts['G']),        
                str(basecnts['T']),
                str(basecnts['N']),
                str(float(ref_cnt)/float(tot_cnt)),
                str(tot_cnt),
                str(pbinom),
                basecnts_hap1['warning'],
                basecnts_hap2['warning']
            ])+'\n')

rmvdhets_file.close()
