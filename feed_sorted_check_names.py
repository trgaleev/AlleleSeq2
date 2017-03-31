import sys
from collections import defaultdict

cnames1   = ['chr' + str(c) for c in range(1,23) + ['X', 'Y', 'M']]
cnames2   = [str(c) for c in range(1,23) + ['X', 'Y', 'MT']]
cnamesbam = []

with open(sys.argv[1],'r') as cf:
   for line in cf:
       cname = line.split()[0]
       if cname in cnames1 or cname in cnames2: 
           cnamesbam.append(cname)
if len(cnamesbam) != 25: sys.exit(sys.argv[0]+': unexpected? chromosome names in '+ sys.argv[2])
      

hets_dict = defaultdict(list)

for line in sys.stdin:
   chrm, c2, c3, c4 = line.split('\t')

   # now deal with chr names; just two common conventions for now
   if chrm not in cnamesbam:

       if   chrm in cnames2 and cnamesbam[0] in cnames1: chrm = 'chr' + chrm
       elif chrm in cnames1 and cnamesbam[0] in cnames2: chrm = chrm[3:]
       elif 'M' in chrm:
           for cname in cnamesbam:
               if 'M' in cname: chrm = cname
       else: sys.exit(sys.argv[0]+': unexpected? chromosome names in stdin') 

   hets_dict[chrm].append('\t'.join([chrm, c2, c3, c4]))   

for c in cnamesbam:
	for l in hets_dict[c]: sys.stdout.write(l)

