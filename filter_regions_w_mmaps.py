import sys
import binom


def outwrite(l, log_mm='.'):
    sys.stdout.write('\t'.join([l.strip(), log_mm])+'\n')
                                   
def logwrite(region, mm_h1_count, mm_h2_count, comment):
    log.write('\t'.join([region, mm_h1_count + '_', comment])+'\n')

def rmsiteswrite(region, mm_h1_count, mm_h2_count, comment):
    rm_regions_f.write('\t'.join(region + [mm_h1_count+'_'+mm_h2_count+'__'+comment])+'\n')

log = open(sys.argv[3], 'w')
rm_regions_f = open(sys.argv[2], 'w')
mode = sys.argv[1]
log.write('\t'.join(['region', 'mm_count:h1_h2', 'mm_log\n']))

mm_counts_dict = {}
with open(sys.argv[4], 'r') as mm_file:
    mm_file.readline()
    for line in mm_file: 
        region, h1_count, h2_count = line.strip().split('\t')
        mm_counts_dict[region] = {'h1_count':h1_count, 'h2_count':h2_count}


sys.stdout.write(sys.stdin.readline().strip()+'\tmmap_log\n')


for line in sys.stdin:
    region, h1_count, h2_count, h1_allele_ratio, p_binom, snv_count = line.strip().split('\t')

    if region not in mm_counts_dict:
        outwrite(l=line, log_mm='no_mmap_reads')
        logwrite(region, '0', '0', 'no_mmap_reads')
        continue

    # now, only bother if the unique counts aren't equal and the 'weaker' allele is seen in multi-mapped reads

    if   int(h1_count) > int(h2_count): 
        if region in mm_counts_dict: weaker_hap_mm_count = int(mm_counts_dict[region]['h2_count'])
        else: weaker_hap_mm_count = 0
        weaker_hap_un_count = int(h2_count)
        weaker_hap = 'h2'
    elif int(h1_count) < int(h2_count): 
        if region in mm_counts_dict: weaker_hap_mm_count = int(mm_counts_dict[region]['h1_count'])
        else: weaker_hap_mm_count = 0
        weaker_hap_un_count = int(h1_count)
        weaker_hap = 'h1'
    else: 
        outwrite(l=line)
        logwrite(region, mm_counts_dict[region]['h1_count'], mm_counts_dict[region]['h2_count'], 'eq_uniq_cnts')
        continue



    if  weaker_hap_mm_count == 0: 
        outwrite(l=line, log_mm=weaker_hap+':0')
        logwrite(line, mm_counts_dict[region]['h1_count'], mm_counts_dict[region]['h2_count'],  weaker_hap+':0')
    else: 
        # now, either remove the site if the allelic ratio changes by more than a threshold
        if mode != 'adjust':
            new =  float(weaker_hap_un_count + weaker_hap_mm_count) / (float(h1_count) + float(h2_count) + float(weaker_hap_mm_count)) 
            old =  float(weaker_hap_un_count) / (float(h1_count) + float(h2_count))
            if (new - old) > float(mode):
                logwrite(region, mm_counts_dict[region]['h1_count'],  mm_counts_dict[region]['h2_count'], '_'.join([str(new-old), weaker_hap, 'removed']))
                rmsiteswrite(region, mm_counts_dict[region]['h1_count'],  mm_counts_dict[region]['h2_count'], weaker_hap+';'+str(new-old))
            else:
                logwrite(region, mm_counts_dict[region]['h1_count'], mm_counts_dict[region]['h2_count'], '_'.join([str(new-old), weaker_hap, 'within_thresh']))
                outwrite(region, weaker_hap + ':' +  str(weaker_hap_mm_count) +';within_thresh')
                

        # or adjust counts in the most conservative way: add the mm counts to the weaker allele:
        # all of them or until balanced with the stonger
        # to make sure the count imbalance is not caused by the multimapping reads              
        else:
            adj = min((weaker_hap_un_count + weaker_hap_mm_count), max(int(h1_count), int(h2_count))) 
            diff = adj - weaker_hap_un_count 
            if weaker_hap == 'h1': h1_count = adj
            else: h2_count = adj
            
            new_tot = int(h1_count) + int(h2_count) 
            new_h1_allele_ratio = float(h1_count)/float(new_tot)
            new_p_binom = binom.binomtest(int(h1_count), new_tot, 0.5)

            sys.stdout.write('\t'.join([
                region, 
                str(h1_count), str(h2_count), 
                str(new_h1_allele_ratio),
                str(new_p_binom), 
                snv_count,
                weaker_hap + ':+' + str(int(diff))
                ])+'\n')

            logwrite(line, mm_counts_dict[region]['h1_count'], mm_counts_dict[region]['h2_count'], weaker_hap+':+'+str(int(diff)))


rm_regions_f.close()
log.close()
