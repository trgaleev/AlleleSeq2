import sys
import binom


def outwrite(l, log_mm='.'):
    sys.stdout.write('\t'.join([l.strip(), log_mm])+'\n')
                                   
def logwrite(region, mm_hap1_count, mm_hap2_count, comment):
    log.write('\t'.join([region, mm_hap1_count + '_', comment])+'\n')

def rmsiteswrite(region, mm_hap1_count, mm_hap2_count, comment):
    rm_regions_f.write('\t'.join(region + [mm_hap1_count+'_'+mm_hap2_count+'__'+comment])+'\n')

log = open(sys.argv[3], 'w')
rm_regions_f = open(sys.argv[2], 'w')
mode = sys.argv[1]
log.write('\t'.join(['region', 'mm_count:hap1_hap2', 'mm_log\n']))

mm_counts_dict = {}
with open(sys.argv[4], 'r') as mm_file:
    mm_file.readline()
    for line in mm_file: 
        region, hap1_count, hap2_count = line.strip().split('\t')
        mm_counts_dict[region] = {'hap1_count':hap1_count, 'hap2_count':hap2_count}


sys.stdout.write(sys.stdin.readline().strip()+'\tmmap_log\n')


for line in sys.stdin:
    region, hap1_count, hap2_count, hap1_allele_ratio, p_binom, snv_count = line.strip().split('\t')

    if region not in mm_counts_dict:
        outwrite(l=line, log_mm='no_mmap_reads')
        logwrite(region, '0', '0', 'no_mmap_reads')
        continue

    # now, only bother if the unique counts aren't equal and the 'weaker' allele is seen in multi-mapped reads

    if   int(hap1_count) > int(hap2_count): 
        if region in mm_counts_dict: weaker_hap_mm_count = int(mm_counts_dict[region]['hap2_count'])
        else: weaker_hap_mm_count = 0
        weaker_hap_un_count = int(hap2_count)
        weaker_hap = 'hap2'
    elif int(hap1_count) < int(hap2_count): 
        if region in mm_counts_dict: weaker_hap_mm_count = int(mm_counts_dict[region]['hap1_count'])
        else: weaker_hap_mm_count = 0
        weaker_hap_un_count = int(hap1_count)
        weaker_hap = 'hap1'
    else: 
        outwrite(l=line)
        logwrite(region, mm_counts_dict[region]['hap1_count'], mm_counts_dict[region]['hap2_count'], 'eq_uniq_cnts')
        continue



    if  weaker_hap_mm_count == 0: 
        outwrite(l=line, log_mm=weaker_hap+':0')
        logwrite(line, mm_counts_dict[region]['hap1_count'], mm_counts_dict[region]['hap2_count'],  weaker_hap+':0')
    else: 
        # now, either remove the site if the allelic ratio changes by more than a threshold
        if mode != 'adjust':
            new =  float(weaker_hap_un_count + weaker_hap_mm_count) / (float(hap1_count) + float(hap2_count) + float(weaker_hap_mm_count)) 
            old =  float(weaker_hap_un_count) / (float(hap1_count) + float(hap2_count))
            if (new - old) > float(mode):
                logwrite(region, mm_counts_dict[region]['hap1_count'],  mm_counts_dict[region]['hap2_count'], '_'.join([str(new-old), weaker_hap, 'removed']))
                rmsiteswrite(region, mm_counts_dict[region]['hap1_count'],  mm_counts_dict[region]['hap2_count'], weaker_hap+';'+str(new-old))
            else:
                logwrite(region, mm_counts_dict[region]['hap1_count'], mm_counts_dict[region]['hap2_count'], '_'.join([str(new-old), weaker_hap, 'within_thresh']))
                outwrite(region, weaker_hap + ':' +  str(weaker_hap_mm_count) +';within_thresh')
                

        # or adjust counts in the most conservative way: add the mm counts to the weaker allele:
        # all of them or until balanced with the stonger
        # to make sure the count imbalance is not caused by the multimapping reads              
        else:
            adj = min((weaker_hap_un_count + weaker_hap_mm_count), max(int(hap1_count), int(hap2_count))) 
            diff = adj - weaker_hap_un_count 
            if weaker_hap == 'hap1': hap1_count = adj
            else: hap2_count = adj
            
            new_tot = int(hap1_count) + int(hap2_count) 
            new_hap1_allele_ratio = float(hap1_count)/float(new_tot)
            new_p_binom = binom.binomtest(int(hap1_count), new_tot, 0.5)

            sys.stdout.write('\t'.join([
                region, 
                str(hap1_count), str(hap2_count), 
                str(new_hap1_allele_ratio),
                str(new_p_binom), 
                snv_count,
                weaker_hap + ':+' + str(int(diff))
                ])+'\n')

            logwrite(line, mm_counts_dict[region]['hap1_count'], mm_counts_dict[region]['hap2_count'], weaker_hap+':+'+str(int(diff)))


rm_regions_f.close()
log.close()
