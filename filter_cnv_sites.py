import sys

pass_dict={}
with open(sys.argv[1], 'r') as cnvfile:
    for line in cnvfile:
        if not line.startswith('chrm'):
            chrm,snppos,rd = line.strip().split('\t')
            if float(rd) <= 1.5 and float(rd) >= 0.5:
                pass_dict[chrm+'_'+snppos] = rd


with open(sys.argv[2], 'r') as countsfile:
    for line in countsfile:
        if not line.startswith('#'):
            if line.split('\t')[0]+"_"+line.split('\t')[1] in pass_dict:
                sys.stdout.write(line.strip()+'\t'+pass_dict[line.split('\t')[0]+"_"+line.split('\t')[1]]+'\n')
        else: sys.stdout.write(line.strip()+'\tcnv\n')
