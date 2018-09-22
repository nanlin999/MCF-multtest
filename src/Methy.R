library(qvalue)  
library(Rcpp)
sourceCpp("~/MCF.cpp")

### New comparison using data with 1 replicate ###

Brain_CpG = read.table("~/Data/Brain/HuFGM02_BrainGerminalMatrix_Bisulfite-Seq_A04698_CpG.bedGraph", header = F)
Brain_density = read.table("~/Data/Brain/HuFGM02_BrainGerminalMatrix_Bisulfite-Seq_A04698_density.bedGraph", header = F)

colnames(Brain_CpG) = c('chr', 'start', 'end', 'ratio')
colnames(Brain_density) = c('chr', 'start', 'end', 'total')

brain_all = merge(Brain_CpG, Brain_density, by = c('chr','start','end'), sort = F)
brain_all[,'brain_methy'] = round(brain_all$ratio * brain_all$total)
brain_all[,'brain_unmethy'] = brain_all$total - brain_all$brain_methy

H1ES = read.table("~/Data/H1ES/E003.methcount", header = F)
colnames(H1ES) = c('chr', 'start', 'end', 'es_methy', 'es_total')
H1ES[,'es_unmethy'] = H1ES$es_total - H1ES$es_methy

skin_cpg = read.table("~/Data/Skin/E058.fm.bedGraph", header = F)
skin_density = read.table("~/Data/Skin/E058.rc.bedGraph", header = F)

colnames(skin_cpg) = c('chr', 'start', 'end', 'ratio')
colnames(skin_density) = c('chr', 'start', 'end', 'total')
skin_all = merge(skin_cpg, skin_density, by = c('chr','start','end'), sort = F)
skin_all[,'skin_methy'] = round(skin_all$ratio * skin_all$total)
skin_all[,'skin_unmethy'] = skin_all$total - skin_all$skin_methy



sub1 = brain_all[,c('chr', 'start', 'end', 'brain_methy', 'brain_unmethy')]
sub2 = skin_all[,c('chr', 'start', 'end', 'skin_methy', 'skin_unmethy')]
sub2 = H1ES[,c('chr', 'start', 'end', 'es_methy', 'es_unmethy')]

Brain_ES_raw = merge(sub1, sub2, by = c('chr','start','end'), sort = F)
Skin_ES_raw = merge(sub1, sub2, by = c('chr','start','end'), sort = F)
Brain_Skin_raw = merge(sub1, sub2, by = c('chr','start','end'), sort = F)


ORG = Brain_Skin_raw
Z1 = Brain_Skin_raw$brain_methy
Z2 = Brain_Skin_raw$brain_unmethy
Z3 = Brain_Skin_raw$skin_methy
Z4 = Brain_Skin_raw$skin_unmethy


# remove small counts
counts_cutoff = 15
idx = which((Z1+Z2)<counts_cutoff|(Z3+Z4)<counts_cutoff)
V = ORG[-idx,]
Z1 = Z1[-idx]
Z2 = Z2[-idx]
Z3 = Z3[-idx]
Z4 = Z4[-idx]

n.test = nrow(V)
p.org <- rep(0,n.test)
p.next <- rep(0,n.test)

for(i in 1:n.test){
  n1 <- Z1[i] + Z2[i]
  n2 <- Z3[i] + Z4[i]
  a <- Z1[i]
  b <- Z3[i]  
  l <- max(0,a+b-n2)
  u <- min(a+b,n1)
  prob.all <- dhyper(l:u,n1,n2,(a+b))
  prob.obs <- dhyper(a,n1,n2,(a+b))
  p.org[i] <- sum(prob.all[which(prob.all<=prob.obs)])
  p.next[i] <- sum(prob.all[which(prob.all<prob.obs)])
  print(i)
}
# rounding problems
p.org[which(p.org>1)] = 1
p.next[which(p.next>1)] = 1


MT <- qvalue(p.org) ## 5 mins
pi0 = MT$pi0 # 1
QVALUE <- MT$qvalues
# sum(QVALUE<0.1)


pi1 <- 1-pi0

n.rep <- 10
n.ecdf <- n.rep*n.test
randp.ecdf <- rep(0,n.ecdf)
for(i in 0:(n.rep-1)){
  randp.ecdf[(i*n.test+1):((i+1)*n.test)] <- runif(n.test,p.next,p.org)
  print(i)  
}

sourceCpp("~/MCF.cpp")
### under FDR nominal level 0.1 ########################
a = 0.00001
b = 0.1
l = 0
R = 0
pFDR = 0
leng = length(randp.ecdf)
while( abs(pFDR-0.1)>0.0000001 ){
  l = (a+b)/2
  cdf <- sum(randp.ecdf<l)/leng
  pFDR = (1-pi1)*l/cdf  
  pFDR
  if(pFDR > 0.1){
    b=l
  }
  if(pFDR < 0.1){
    a=l
  }  
  print(c(a,b,l,pFDR))
}
print(c(a,b,l,pFDR))
mcf <- MCF(p.org, p.next, l)
rej.0.1 <- which(mcf >= quantile(mcf, probs = 1-cdf, type = 1))
# rej.0.1 <- which(mcf>sort(mcf)[n.test-R])
rejC.0.1 <- which(QVALUE<0.1)

### under FDR nominal level 0.05 ###################
a = 0.00001
b = 0.1
l = 0
R = 0
pFDR = 0
leng = length(randp.ecdf)
while( abs(pFDR-0.05)>0.0000001 ){
  l = (a+b)/2
  cdf <- sum(randp.ecdf<l)/leng
  pFDR = (1-pi1)*l/cdf  
  R = round(n.test*cdf) 
  pFDR
  if(pFDR > 0.05){
    b=l
  }
  if(pFDR < 0.05){
    a=l
  }  
  print(c(a,b,l,pFDR))
}
print(c(a,b,l,pFDR))
mcf <- MCF(p.org, p.next, l)
rej.0.05 <- which(mcf >= quantile(mcf, probs = 1-cdf, type = 1))
rejC.0.05 <- which(QVALUE<0.05)


REJ010 <- rep(0,n.test)
REJC010 <- rep(0,n.test)
REJ010[rej.0.1] <- 1
REJC010[rejC.0.1] <- 1
REJ005 <- rep(0,n.test)
REJC005 <- rep(0,n.test)
REJ005[rej.0.05] <- 1
REJC005[rejC.0.05] <- 1

c(sum(REJ010), sum(REJC010), sum(REJ010*REJC010))
c(sum(REJ005), sum(REJC005), sum(REJ005*REJC005))





ALL = data.frame(cbind(V,p.org,p.next,REJ010,REJC010,REJ005,REJC005))
# colnames(ALL) = c('chrome','start','end','A1','A2','A3','A4','B1','B2','B3','B4', 'porg','pnext', 'MCF010','Qvalue010', 'MCF005','Qvalue005')
colnames(ALL) = c('chrome','start','end','A1','A2','B1','B2', 'porg','pnext', 'MCF010','Qvalue010', 'MCF005','Qvalue005')
write.table(ALL,"~/Skin_ES/ALL.txt",row.name=F,quote=FALSE)

