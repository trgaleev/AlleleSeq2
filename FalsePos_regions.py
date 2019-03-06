# as in AlleleSeq_v1.2a; minor edits

''' Some notes on what is going on here.
Basically, we want to use simulation to explicitly calculate a FDR for binomial tests on unbalanced alleles.  We use
a binomial pvalue test to determine whether the ratio of alleles departs significantly from what would be expected
from a fair coin toss.

However, since the number of trials varies a lot from one test to another, it seems best to use an explicit method.

Imagine that you want a particular overall FDR, say 0.1, for the entire dataset.  Then, what pvalue threshhold would correspond to that?

say we have n trials, where the number of coin flips in a trial varies, and is given by cnt(i)

FDR = Nsim/Nact, where:

Nsim = sum( indicator(test(i) < pval)) over i.  This is the number of trials of the fair coin that had a "surprising" outcome, i.e. 
were further in the tail than the pval threshold.  In a perfect, non-discrete world, Nsim/n would equal pval, but the whole point of this
exercise is that in the discrete, imperfect world it doesnt.

Nact = the number of actual datapoints observed to have a binomial probability less than pval.

So, for a given pval, say 0.05, we can calculate the FDR, which will be larger.  The first output from this code consists of a nice sampling of
example pvals and their corresponding FDR.  We are interested in the reverse of this, i.e. having picked an FDR, we want the pval that would best give us
this FDR.

Thats the point of the second part of the output.  Starting from the largest pval, we work our way down, calculating FDR for each test,
until FDR falls below our target.

Note that FDR is NOT monotonically decreasing as we do this. Its true that both Nsim and Nact decrease.  However, Nact is strictly decreasing, but Nsim can hold steady, which results in temporarily increasing FDR over that interval.

Also note that we do multiple simulations and use the mean of the Nsim values, in order to smooth out the results.

'''
## import sys, getInterestingHetsAnnotations, bisect, binom, random, numpy, pdb
import sys, bisect, binom, random, numpy, pdb

class binomMemo(object):
    def __init__(self, n):
        self.n=n
        self.cache=[[binom.binomtest(j, i, 0.5) for j in range(i+1)] for i in range(n)]
    def binomtest(self, a, cnt):
        if cnt<self.n:
            return self.cache[cnt][a]
        else:
            return binom.binomtest(a, cnt, 0.5)

def simpval(cnt,bm):
    a=sum([random.randint(0,1) for i in range(cnt)])
    pval=bm.binomtest(a, cnt)
    return pval

def simpval2(cnt,bm):
    a=sum([random.randint(0,1) for i in range(cnt)])
    pval=bm.binomtest(a, cnt)
    return a

if __name__=='__main__':
    ifile=sys.argv[1]
    sims=int(sys.argv[2])
    #verbose=False
    verbose = len(sys.argv)==5 and sys.argv[4]=='-v'
    bestFDR=bestPV=None

    random.seed(0)
    target=float(sys.argv[3]) # target is the FDR we are looking for, we want to find the corresponding pval

    print "#"," ".join(sys.argv)
    print "pval\tP\tFP\tFDR"
    bm=binomMemo(60)

## let's try to make it simpler, without the additional module
##    h=getInterestingHetsAnnotations.Handler(ifile, hasHeader=True)
##    n=h.getCount()
##    g=h.getAllAnnotationsGenerator();
##    act_pvals=numpy.zeros(n) # pval as reported in counts file
##    cnt_sums=numpy.zeros(n, dtype=numpy.int)  # sum of major and minor alleles
##
##    # for each hetSNP, get the count of major and minor allele from the input file
##    for i, t in enumerate(g):
##        c,pos,rec = t
##        d=h.toDict(rec)
##        act_pvals[i]=float(d['SymPval'])
##        counts=[d['cA'],d['cC'],d['cG'],d['cT']]
##        counts=[int(e) for e in counts]
##        counts=sorted(counts, reverse=True)[0:2]
##        cnt_sums[i]=sum(counts)

    cnt_sums_list  = []
    act_pvals_list = []
    with open(ifile, 'r') as inf:
	inf.readline()
	for line in inf:
            region,hap1_count,hap2_count,hap1_allele_ratio,p_binom,snv_count,snv_coords,mmap_log = line.split('\t')
            act_pvals_list.append(float(p_binom))
            #counts = [int(e) for e in [cA, cC, cG, cT]]
            #counts = sorted(counts, reverse=True)[:2]
            counts = [int(hap1_count), int(hap2_count)]
            cnt_sums_list.append(sum(counts))
    act_pvals = numpy.array(act_pvals_list)
    cnt_sums  = numpy.array(cnt_sums_list)
    n = len(cnt_sums)

##
      
	
            	
    
    act_pvals.sort()
    

    sim_pvals=numpy.array([ sorted([simpval(cnt_sums[j],bm) for j in xrange(n)]) for i in xrange(sims)])
   
    #sim_pvals_means=numpy.mean(sim_pvals, 0)

    pvs=[e*0.001 for e in range(10)]+[e*0.01 for e in range(1,10)]+[e*0.1 for e in range(1,10)]
    # for a given test pv, find the number of actual pvals that are smaller, and the number of sim pvals that are smaller.
    # FDR is the ratio
    for pv in pvs:
        Nact=bisect.bisect(act_pvals, pv)
        mean_Nsims=numpy.mean([bisect.bisect(sim_pvals[i], pv) for i in xrange(sims)])
        FDR=mean_Nsims/(Nact+1)
        print "%f\t%s\t%f\t%f" % (pv, Nact, mean_Nsims, FDR)

    # This is my attempt to find the act_pval that corresponds best to the desired target FDR.  
    # This version walks from largest observed pvalue to the smallest.
    if target:
        last_FDR=last_pv=0.0
        for Nact, pv in sorted(enumerate(act_pvals), reverse=True):
            mean_Nsims=numpy.mean([bisect.bisect(sim_pvals[i], pv) for i in xrange(sims)])
            FDR=mean_Nsims/(Nact+1)
            if verbose: print "test %d %f %f %f" % (Nact,mean_Nsims,FDR, pv)
#tg20161115: this seems to fail (in extreme cases) if the first bestFDR < target is 0
#            if not bestFDR and FDR < target:
            if bestFDR==None and FDR < target:
                print "target %f" % target
                print "before %f %f" % (last_FDR, last_pv)
                print "after  %f %f" % (FDR, pv)
                bestFDR = FDR; bestPV = pv
#tg20161115: stop tests if bestFDR found, bestFDR and bestPV don't change anyway:
		if not verbose: break

            last_FDR=FDR; last_pv=pv

#tg20161115: to avoid crashing when bestFDR and bestPV are None. 
	#print "Target %f FDR %f pv %f" % (target,bestFDR, bestPV)
	if bestFDR!=None and bestPV!=None: print "Target %f FDR %f pv %f" % (target,bestFDR, bestPV)
	else: 
		print "# couldn't calculate pvalue threshold: no snvs in input or none with small enough p-value (?)\n# will set p value threshold to 0: \nTarget "+str(target)+" FDR NA pv 0" 
		print >> sys.stderr, "WARNING "+sys.argv[0]+": couldn't calculate pvalue threshold: no snvs in input or none with small enough p-value (?). will set p value threshold to 0 - no snvs will be reported after filtering" 

