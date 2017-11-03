import sys
import read_pileup
import binom

def outwrite(l, log_mm='.'):
    sys.stdout.write('\t'.join([l.strip(), log_mm])+'\n')
                                   
def logwrite(l, mm_counts, comment):
    mm_cnts='_'.join([str(x) for x in (mm_counts['A'],mm_counts['C'],mm_counts['G'],mm_counts['T'])])
    log.write('\t'.join(l.split('\t')[:2] + [mm_cnts, comment])+'\n')

def rmsiteswrite(l, mm_counts, comment):
    mm_cnts='_'.join([str(x) for x in (mm_counts['A'],mm_counts['C'],mm_counts['G'],mm_counts['T'])])
    rm_hetSNV_f.write('\t'.join(l.split('\t')[:2] + [mm_cnts+'__'+comment])+'\n')


log = open(sys.argv[3],'w')
rm_hetSNV_f = open(sys.argv[2], 'a')
mode = sys.argv[1]
log.write('\t'.join(['chr', 'pos', 'max_mm_count:A_C_G_T_N', 'mm_h1_warn;mm_h2_warn;mm_log\n']))

mm_pileup_dict = read_pileup.pileup_to_basecnts(sys.argv[4:])

sys.stdout.write(sys.stdin.readline().strip()+'\tmmap_log\n')

for line in sys.stdin:
    (chrm, ref_coord, h1_coord, h2_coord, ref_allele, h1_allele, h2_allele, 
            cA, cC, cG, cT, cN, ref_allele_ratio, sum_ref_n_alt_cnts, 
            p_binom, cnv) = line.split('\t')
  
    h1_mm_basecnts = mm_pileup_dict.get(h1_coord,{'A':0, 'C':0, 'G':0, 'T':0, 'N':0, 'warning':'.'})
    h2_mm_basecnts = mm_pileup_dict.get(h2_coord,{'A':0, 'C':0, 'G':0, 'T':0, 'N':0, 'warning':'.'})

    # again because of the possibility of reads flanking, say, a misphased indel, will check *both* alleles for multi-mapping reads bearing each allele
    mm_basecnts = {
            'A': max(h1_mm_basecnts['A'], h2_mm_basecnts['A']),
            'C': max(h1_mm_basecnts['C'], h2_mm_basecnts['C']),
            'G': max(h1_mm_basecnts['G'], h2_mm_basecnts['G']),
            'T': max(h1_mm_basecnts['T'], h2_mm_basecnts['T']),
            'N': max(h1_mm_basecnts['N'], h2_mm_basecnts['N'])
    }

    # now, only bother if the unique counts aren't equal and the 'weaker' allele is seen in multi-mapped reads
    um_basecnts = {'A':float(cA), 'C':float(cC), 'G':float(cG), 'T':float(cT), 'N':float(cN)}
    if   um_basecnts[h1_allele] > um_basecnts[h2_allele]: 
        weaker = h2_allele
    elif um_basecnts[h1_allele] < um_basecnts[h2_allele]: 
        weaker = h1_allele
    else: 
        outwrite(l=line)
        logwrite(line, mm_basecnts, ';'.join([h1_mm_basecnts['warning'], h2_mm_basecnts['warning'], 'eq_uniq_cnts']))
        continue



    if mm_basecnts[weaker] == 0: 
        outwrite(l=line, log_mm=weaker+':'+str(mm_basecnts[weaker]))
        logwrite(line, mm_basecnts, ';'.join([h1_mm_basecnts['warning'], h2_mm_basecnts['warning'], weaker+':'+str(mm_basecnts[weaker])]))
    else: 
        # now, either remove the site if the allelic ratio changes by more than a threshold
        if mode != 'adjust':
            new = (um_basecnts[weaker] + float(mm_basecnts[weaker])) / (um_basecnts[h1_allele] + um_basecnts[h2_allele] + mm_basecnts[weaker]) 
            old = um_basecnts[weaker] / (um_basecnts[h1_allele] + um_basecnts[h2_allele])
            if (new - old) > float(mode):
               #logwrite(line, mm_basecnts, ';'.join([h1_mm_basecnts['warning'], h2_mm_basecnts['warning'],'.', weaker+'_removed']))
                logwrite(line, mm_basecnts, ';'.join([h1_mm_basecnts['warning'], h2_mm_basecnts['warning'], weaker+'_removed']))
                rmsiteswrite(line, mm_basecnts, weaker+';')
            else:
               #logwrite(line, mm_basecnts, ';'.join([h1_mm_basecnts['warning'], h2_mm_basecnts['warning'],'.', weaker+'_within_thresh']))
                logwrite(line, mm_basecnts, ';'.join([h1_mm_basecnts['warning'], h2_mm_basecnts['warning'], weaker+'_within_thresh']))
                outwrite(line, weaker+':'+str(mm_basecnts[weaker])+';within_thresh')
                

        # or adjust counts in the most conservative way: add the mm counts to the weaker allele:
        # all of them or until balanced with the stonger
        # to make sure the count imbalance is not caused by the multimapping reads              
        else:
            adj = min((um_basecnts[weaker] + mm_basecnts[weaker]), max(um_basecnts[h1_allele], um_basecnts[h2_allele])) 
            diff = adj - um_basecnts[weaker] 
            um_basecnts[weaker] = adj
            
            new_tot = int(um_basecnts[h1_allele] + um_basecnts[h2_allele])
            new_ratio = float(um_basecnts[ref_allele])/float(new_tot)
            new_p_binom = binom.binomtest(um_basecnts[ref_allele], new_tot, 0.5)

            sys.stdout.write('\t'.join([
                chrm, ref_coord, h1_coord, h2_coord, ref_allele, h1_allele, h2_allele,
                str(int(um_basecnts['A'])),
                str(int(um_basecnts['C'])),
                str(int(um_basecnts['G'])),
                str(int(um_basecnts['T'])),
                str(int(um_basecnts['N'])),
                str(new_ratio), str(new_tot), str(new_p_binom), 
                cnv[:-1], weaker+':+'+str(int(diff))
                ])+'\n')

            logwrite(line, mm_basecnts, ';'.join([h1_mm_basecnts['warning'], h2_mm_basecnts['warning'], weaker+':+'+str(int(diff))]))


rm_hetSNV_f.close()
log.close()
