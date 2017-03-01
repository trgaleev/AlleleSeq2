"""
http://www.htslib.org/doc/samtools.html

"In the pileup format (without -u or -g), each line represents a genomic position, consisting of chromosome name, 1-based coordinate, reference base, the number of reads covering the site, read bases, base qualities and alignment mapping qualities. Information on match, mismatch, indel, strand, mapping quality and start and end of a read are all encoded at the read base column.

At this column, a dot stands for a match to the reference base on the forward strand, a comma for a match on the reverse strand, a '>' or '<' for a reference skip, `ACGTN' for a mismatch on the forward strand and `acgtn' for a mismatch on the reverse strand. A pattern `\\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion between this reference position and the next reference position. The length of the insertion is given by the integer in the pattern, followed by the inserted sequence. Similarly, a pattern `-[0-9]+[ACGTNacgtn]+' represents a deletion from the reference. The deleted bases will be presented as `*' in the following lines. Also at the read base column, a symbol `^' marks the start of a read. The ASCII of the character following `^' minus 33 gives the mapping quality. A symbol `$' marks the end of a read segment."
"""


import sys

def pileup_to_basecnts (filelist):
    pileup_dict = {}

    for mf in filelist:
        with open(mf,'r') as in_m:
            for line in in_m:
                warning = '.'
                chr, crd, a, tot_pileup_cnt, seq, _ = line.split()
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
                if (basecnts['A'] + basecnts['C'] + basecnts['G'] + basecnts['T'] + basecnts['N'] + deleted_bases + spliced) != int(tot_pileup_cnt):
                    sys.exit('error1: unexpected base counts / symbols in pileup line\n'+line)

                if sum(basecnts.values()) > basecnts[a]: warning = str(sum(basecnts.values()) - basecnts[a])+'_other_alleles'                
                if max(basecnts.values()) > basecnts[a]: warning = 'hap_allele_not_largest_cnt'
                if sum(basecnts.values()) == 0: warning = 'zero_cnt'

                basecnts['warning'] = warning

                pileup_dict[chr+'_'+crd] = basecnts

    return pileup_dict
