"""
http://www.htslib.org/doc/samtools.html

In the pileup format (without -u or -g), each line represents a genomic position, consisting of chromosome name, 1-based coordinate, reference base, the number of reads covering the site, read bases, base qualities and alignment mapping qualities. Information on match, mismatch, indel, strand, mapping quality and start and end of a read are all encoded at the read base column.

At this column, a dot stands for a match to the reference base on the forward strand, a comma for a match on the reverse strand, a '>' or '<' for a reference skip, `ACGTN' for a mismatch on the forward strand and `acgtn' for a mismatch on the reverse strand. A pattern `\\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion between this reference position and the next reference position. The length of the insertion is given by the integer in the pattern, followed by the inserted sequence. Similarly, a pattern `-[0-9]+[ACGTNacgtn]+' represents a deletion from the reference. The deleted bases will be presented as `*' in the following lines. Also at the read base column, a symbol `^' marks the start of a read. The ASCII of the character following `^' minus 33 gives the mapping quality. A symbol `$' marks the end of a read segment.
"""


import binom
import string
import os
import sys
import scipy.stats
import re
from collections import defaultdict


hetSNV_dict=defaultdict(dict)
hetSNV_list=[]

# read the h1 and h2 coords
with open(sys.argv[1],'r') as in_h2:
    for line in in_h2:
        line = line.strip()
        h_chr, _, h_crd, ref_info = line.split('\t')
        #h_cr, h = h_chr.split('_')
        r_chr, r_crd, r_a, a1, a2 = ref_info.split('_')
        #hetSNV_dict[r_chr+'_'+r_crd][h]={'crd':h_crd}
        hetSNV_dict[r_chr+'_'+r_crd]['h1_pos'] = h_chr+'_'+h_crd
        hetSNV_dict[r_chr+'_'+r_crd]['r_a'] = r_a
        hetSNV_dict[r_chr+'_'+r_crd]['a1'] = a1
        hetSNV_dict[r_chr+'_'+r_crd]['a2'] = a2

with open(sys.argv[2],'r') as in_h2:
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
mpileup_dict = {}
mpileup_comment_dict = {}
for mf in sys.argv[3:5]:
    with open(mf,'r') as in_m:
        for line in in_m:
            comment = '.'
            chr, crd, a, tot_mpileup_cnt, seq, _ = line.split()
            a = a.upper()
            basecnts={'A':0, 'C':0, 'G':0, 'T':0, 'N':0}

            # remove indications of insertions and deletions b/w current and next pos
            if '+' in seq or '-' in seq:
                tmp_seq = ''
                new_indel = False
                number = ''
                for character in seq:
                    if not new_indel:
                        if character not in ['+','-']: tmp_seq += character
                        else: 
                            new_indel = True
                            number = ''
                    else:
                        if character.isalpha(): 
                            number = int(number)
                            number -= 1
                        else: number += character
                        if int(number) == 0: new_indel = False
                    
                seq = tmp_seq



            basecnts[a] = seq.count(',') + seq.count('.')
            deleted_bases = seq.count('*')
            spliced = seq.count('<') + seq.count('>')

            for base in basecnts: 
                if base in seq.upper(): basecnts[base] = seq.upper().count(base)
            if (basecnts['A'] + basecnts['C'] + basecnts['G'] + basecnts['T'] + basecnts['N'] + deleted_bases + spliced) != int(tot_mpileup_cnt): 
                sys.exit('error1: unexpected base counts / symbols in mpileup line\n'+line)
            

            if max(basecnts.values()) > basecnts[a]: comment = 'allele_not_max'

            mpileup_dict[chr+'_'+crd] = basecnts
            mpileup_comment_dict[chr+'_'+crd] = comment

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
    'p_binom',
    'note_h1',
    'note_h2'
])+'\n')

#for k in hetSNV_dict:
# to keep the original order
for k in hetSNV_list:

    basecnts_h1 = mpileup_dict.get(hetSNV_dict[k].get('h1_pos',None), {'A':0, 'C':0, 'G':0, 'T':0, 'N':0})
    basecnts_h2 = mpileup_dict.get(hetSNV_dict[k].get('h2_pos',None), {'A':0, 'C':0, 'G':0, 'T':0, 'N':0})
    comment_h1 = mpileup_comment_dict.get(hetSNV_dict[k].get('h1_pos',None),'.')
    comment_h2 = mpileup_comment_dict.get(hetSNV_dict[k].get('h2_pos',None),'.')
    
    basecnts = {
        'A':basecnts_h1['A']+basecnts_h2['A'],
        'C':basecnts_h1['C']+basecnts_h2['C'],
        'G':basecnts_h1['G']+basecnts_h2['G'],
        'T':basecnts_h1['T']+basecnts_h2['T'],
        'N':basecnts_h1['N']+basecnts_h2['N']
    }

    tot_cnt = basecnts[hetSNV_dict[k]['a1']]+basecnts[hetSNV_dict[k]['a2']]

    if tot_cnt >= int(sys.argv[5]):

        ref_cnt = basecnts[hetSNV_dict[k]['r_a']]
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
            str(pbinom),
            comment_h1,
            comment_h2
        ])+'\n')

