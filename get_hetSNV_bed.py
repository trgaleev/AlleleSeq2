import sys
from collections import defaultdict

cnames1   = ['chr' + str(c) for c in range(1,23) + ['X', 'Y', 'M']]
cnames2   = [str(c) for c in range(1,23) + ['X', 'Y', 'MT']]
cnamesref = []

with open(sys.argv[2],'r') as cf:
   for line in cf:
       if line.startswith('>'):
           cname = line.strip()[1:]
           if cname in cnames1 or cname in cnames2: 
               cnamesref.append(cname)
if len(cnamesref) != 25: sys.exit(sys.argv[0]+': unexpected? chromosome names in '+ sys.argv[2])
      

vcf_hets_dict = defaultdict(list)

for line in sys.stdin:
    if line.startswith('#CHROM'): sample_col_idx = line.split().index(sys.argv[1])
    if not line.startswith('#'):
        CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT = line.split('\t')[:9]
        if FILTER=='PASS': 
            GTl=line.split('\t')[sample_col_idx].split(':')[0].replace('/','|').split('|')
            al=[REF]+ALT.split(',')
            if al[int(GTl[0])]!=al[int(GTl[1])] and len(al[int(GTl[0])])==len(al[int(GTl[1])])==len(REF)==1:

                # now deal with chr names; just two common conventions for now
                if CHROM not in cnamesref:

                    if   CHROM in cnames2 and cnamesref[0] in cnames1: CHROM = 'chr' + CHROM
                    elif CHROM in cnames1 and cnamesref[0] in cnames2: CHROM = CHROM[3:]
                    elif 'M' in CHROM: 
                        for cname in cnamesref:
                            if 'M' in cname: CHROM = cname
                    else: continue # don't write 

                sys.stdout.write('\t'.join([   
                                              CHROM,
                                              str(int(POS)-1),
                                              POS,
                                              '_'.join([CHROM,POS,REF, al[int(GTl[0])], al[int(GTl[1])]])        
                                               ])+'\n')


