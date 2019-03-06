#20170212 cjieming/alleleDB alleledb_v2.0 script
#Chen J et al. (2016). Nat Commun. 7:11101 
# with minor modifications to incorporate into the pipeline



## allow cmd line arguments
args=(commandArgs(TRUE))

#folder=args[1]
if (length(args)>5){
	FDR.thresh = as.numeric(args[6])
}else{
	FDR.thresh = 0.05
}




## start script here
#library(VGAM)
library(VGAM, lib.loc="~/R_libs/")

### data

data1 = read.table(args[1], header=T, stringsAsFactors=F, comment.char="%", check.names=F)
#print(head(data1))
#filename2 = "betabinomial/b_chosen.grad.txt"
filename2 = paste0(args[2],'/',"b_chosen.grad.txt")
data2 = read.table(filename2, header=T, stringsAsFactors=F)

## parameters
p=0.5 ## binomial null p

## finding the second highest value for the SNP
cACGT = data.frame(cA=data1$cA,cC=data1$cC,cG=data1$cG,cT=data1$cT)
lower = apply(cACGT,1,function(x) sort(x, partial=3)[3])

## find total num = highest + second highest
higher = apply(cACGT,1,max)
total = higher+lower


## empirical tests
b = data2$b.choice

## VGAM says that "the value 0 is sometimes supported (if so then it corresponds to the usual binomial distribution)"
## to prevent the code from breaking, a pseudo-zero number is added
## and it's essentially the same as the binomial distribution up to 7 dp
if(b==0){  b=1e-10 }

p.bin = apply(data.frame(2 * mapply(pbinom,lower,total,p)),1,function(x) min(x,1))
p.betabin = apply(data.frame(2 * mapply(pbetabinom,lower,total,p,b)),1,function(x) min(x,1))
data1$p.betabin = p.betabin

## simulations
step = 0.0001
p.thresh = data.frame( c(seq(0,0.01,by=0.001), seq(0.01,0.1,by=0.01)[-1], seq(0.1,1,by=0.1)[-1]) ) #
cutoff <- function(x,y) sum(y<=x)

## calc fp from null and empirical counts
fp <- function(w,p,p.thresh,distrib="binomial",b=0)
{
  ## doing the distribution; as.integer converts table entities to integers
  a=lapply(as.integer(w[,1]),function(x) seq(0,x))
  
  if(distrib == "binomial")
  { 
    b = lapply(a,function(x) apply(as.data.frame(2*pbinom(x,max(x),p)),1,function(x) min(x,1)))
  }
  else if(distrib == "betabinomial")
  { 
    b = lapply(a,function(x) apply(as.data.frame(2*pbetabinom(x,max(x),p,b)),1,function(x) min(x,1)))
  }
  
  ## find which ones are below threshold u
  d = lapply(b,function(x) x<=p.thresh)
  
  ## weight them by actual counts
  e = mapply(function(x,y,z) x*y*z, d, b, w[,2])
  
  ## sum up
  f = sapply(e,max)
  g = sum(f)
  
  return(g)
}

# table of empirical counts
w = as.data.frame(table(total), stringsAsFactors=F)
w = w[as.numeric(w[,1]) >=6,] 

# FP
fp.binomial = apply(p.thresh,1,function(x) fp(w,p,x,"binomial"))
fp.binomial = as.data.frame(cbind(p.thresh,fp.binomial))
colnames(fp.binomial) = c("pval","FP.bin")

fp.betabinomial = apply(p.thresh,1,function(x) fp(w,p,x,"betabinomial",b))
fp.betabinomial = as.data.frame(cbind(p.thresh,fp.betabinomial))
colnames(fp.betabinomial) = c("pval","FP.betabin")

## FDR.txt
tp.bin = apply(p.thresh,1,cutoff,y=p.bin)+1
tp.betabin = apply(p.thresh,1,cutoff,y=p.betabin)+1
fdr.bin = fp.binomial[,2] / tp.bin
fdr.betabin = fp.betabinomial[,2] / tp.betabin
p.choice.bin = max(p.thresh[,1][fdr.bin<=FDR.thresh])
p.choice.betabin = max(p.thresh[,1][fdr.betabin<=FDR.thresh])

fdr.choice.bin = max(fdr.bin[fdr.bin<=FDR.thresh])
fdr.choice.betabin = max(fdr.betabin[fdr.betabin<=FDR.thresh])

## bisection method to find p value
bisect <- function(p,p.sim,p.choice,fdr,fdr.threshold,by,distrib="binomial",b=0,w,p.thresh)
{
  p.fdr.e = matrix(0,100,3)
  e.prev = 10
  flag = 3
  ctr = 1
  p.fdr.e[ctr,1] = p.choice
  p.fdr.e[ctr,2] = fdr
  p.fdr.e[ctr,3] = e.prev
  
  while(flag)
  {
    start = max(0,(p.choice - by/2))   
    end = p.choice + by/2
    by = by/4
    
    if(start==0){ start = 5e-4 } ## do not make it 0
    
    range = seq(start,end,by)
    
    for (i in range)
    {
      tp = cutoff(i,p)
      
      if(distrib == "binomial")
      {
        fp = fp(w,p.thresh,i,"binomial")
      }
      else if(distrib == "betabinomial")
      {
        fp = fp(w,p.thresh,i,"betabinomial",b)
      }
      
      fdr.ind = fp/tp
      e.curr = fdr.threshold - fdr.ind
      ctr = ctr + 1
      
      p.fdr.e[ctr,1] = i
      p.fdr.e[ctr,2] = fdr.ind
      p.fdr.e[ctr,3] = e.curr
      e.prev = p.fdr.e[(ctr-1),3]
      p.choice = i
      
      print(paste("i=",i,"start|end|by",start,end,by))
      print(paste("fdr.threshold=",fdr.threshold,"fdr.ind=",fdr.ind,"e.curr=",e.curr)) ##debug
      print(paste("fp=",fp,"tp=",tp))
      if(e.curr < 0){ break }
      
#        print(paste(start,"|",end,"|",i,"|",ctr,"|",by)) ##debug
#        print(paste("fdr.thresh=",fdr.threshold,"|fdr=",fdr.ind,"|fdrmatrix=",p.fdr.e[(ctr-1),2],
#                    "e.curr=",p.fdr.e[ctr,3],"|e.prev=",p.fdr.e[ctr-1,3])) ##debug
    }
#      print(paste(start,"|",end,"|",i,"|",ctr,"|",by)) ##debug
#      print(paste("fdr.thresh=",fdr.threshold,"|fdr=",p.fdr.e[ctr,2],"|fdrprev=",p.fdr.e[ctr-1,2])) 
#      break##debug
#      print(paste("tp=",tp,"|fp=",fp)) ##debug
#      print(paste("e.curr=",e.curr,"|e.prev=",e.prev)) ##debug
    
    if(signif(p.fdr.e[ctr-1,3],3) == signif(p.fdr.e[ctr,3],3)){ flag = 0 }
  }  
  return(p.fdr.e)
}
p.choice.bin.1 = as.data.frame(bisect(p.bin,fp.bin[,2],p.choice.bin,fdr.choice.bin,FDR.thresh,step,"binomial",b=0,w,p))
p.choice.betabin.1 = as.data.frame(bisect(p.betabin,fp.betabinomial[,2],p.choice.betabin,fdr.choice.betabin,FDR.thresh,step,"betabinomial",b,w,p))

p.choice.bin.1 = p.choice.bin.1[p.choice.bin.1[,3]>0,]
p.choice.bin.2 = p.choice.bin.1[nrow(p.choice.bin.1),1]
p.choice.betabin.1 = p.choice.betabin.1[p.choice.betabin.1[,3]>0,]
p.choice.betabin.2 = p.choice.betabin.1[nrow(p.choice.betabin.1),1]

## formatting FDR.txt
FDR.txt = data.frame(cbind(p.thresh,tp.bin,fp.binomial[,2],fdr.bin,
                           tp.betabin,fp.betabinomial[,2],fdr.betabin))
colnames(FDR.txt) <- c("pval","P.bin","FP.bin","FDR.bin",
                       "P.betabin","FP.betabin","FDR.betabin")
FDR.txt[is.na(FDR.txt)] <- 0
FDR.txt[FDR.txt == "Inf"] <- 0 
FDR.txt = FDR.txt[-nrow(FDR.txt),] 


## take in counts.txt and filter by p.betabin and cnv
#interestingHets.betabinom = data1[data1$p.betabin<=p.choice.betabin,]
interestingHets.betabinom = data1[(data1$p.betabin<=p.choice.betabin.2) & (data1$cnv>=0.5 & data1$cnv<=1.5),]
#print (head(interestingHets.betabinom))

## printing files
write.table(data1,file=args[3], sep="\t",
            row.names=FALSE,quote=FALSE)
write.table(interestingHets.betabinom,file=args[4], sep="\t",
            row.names=FALSE,quote=FALSE)
write.table(FDR.txt,file=args[5],sep="\t",row.names=FALSE,quote=FALSE)
write(paste("FDR.threshold =",FDR.thresh),
      file=args[5],append=TRUE)
write(rbind(paste("p.choice.bin.old =",p.choice.bin),paste("p.choice.betabin.old =",p.choice.betabin)),
      file=args[5],append=TRUE)
write(rbind(paste("p.choice.bin =",p.choice.bin.2),paste("p.choice.betabin =",p.choice.betabin.2)),
      file=args[5],append=TRUE)

