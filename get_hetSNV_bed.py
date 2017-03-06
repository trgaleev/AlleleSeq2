import sys

allowed_chr = [str(c) for c in range(1,23) + ['X', 'Y', 'M', 'MT']] + ['chr' + str(c) for c in range(1,23) + ['X', 'Y', 'M', 'MT']]

with sys.stdin as in_vcf:
    for line in in_vcf:
        if line.startswith('#CHROM'): sample_col_idx = line.split().index(sys.argv[1])
        if not line.startswith('#'):
            CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT=line.split('\t')[:9]
            if FILTER=='PASS' and CHROM in allowed_chr:
                GTl=line.split('\t')[sample_col_idx].split(':')[0].replace('/','|').split('|')
                al=[REF]+ALT.split(',')
                if al[int(GTl[0])]!=al[int(GTl[1])] and len(al[int(GTl[0])])==len(al[int(GTl[1])])==len(REF)==1:
                    sys.stdout.write('\t'.join([    
                                                   CHROM,
                                                   str(int(POS)-1),
                                                   POS,
                                                   '_'.join([CHROM,POS,REF, al[int(GTl[0])], al[int(GTl[1])]])        
                                               ])+'\n')
            
            
