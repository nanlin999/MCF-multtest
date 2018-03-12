 vvlibrary(limma)
library(qvalue)
library(edgeR)

all.cd4_es <- read.table("/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD4_ES/ALL.txt", header=T)
all.cd8_es <- read.table("/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD8_ES/ALL.txt", header=T)

rt.cd4.1 <- log2(all.cd4_es$A1/all.cd4_es$A2)
rt.cd4.2 <- log2(all.cd4_es$A3/all.cd4_es$A4)

rt.cd8.1 <- log2(all.cd8_es$A1/all.cd8_es$A2)
rt.cd8.2 <- log2(all.cd8_es$A3/all.cd8_es$A4)

rt.es.1 <- log2(all.cd4_es$B1/all.cd4_es$B2)
rt.es.2 <- log2(all.cd4_es$B3/all.cd4_es$B4)

dat.mtx = matrix(c(rt.cd4.1, rt.cd4.2, rt.es.1, rt.es.2), ncol = 4)
design <- cbind(Grp1=1,Grp2vs1=c(0,0,1,1))

ft = lmFit(dat.mtx, design)
eb.fit <- eBayes(ft)
limma.pvalue = eb.fit$p.value[,2]
write.csv(pvalue, '/tmp/pvalue')

pvalue = read.csv('/tmp/pvalue')[,2]
limma.fdr = p.adjust(limma.pvalue, method = 'BH')

write.csv(limma.fdr, '/tmp/limma.fdr')

rej.limma = ifelse(limma.fdr < 0.1, 1, 0)
rej.limma[is.na(rej.limma)] = 0

rej.mcf = all.cd4_es$MCF010

sum(rej.limma * rej.mcf)
sum(rej.limma)
sum(rej.mcf)









all.cd4_es <- read.table("/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD4_ES/ALL.txt", header=T)

me.cd4.1 <- all.cd4_es$A1
un.cd4.1 <- all.cd4_es$A2 - all.cd4_es$A1

me.cd4.2 <- all.cd4_es$A3
un.cd4.2 <- all.cd4_es$A4 - all.cd4_es$A3

me.es.1 <- all.cd4_es$B1
un.es.1 <- all.cd4_es$B2 - all.cd4_es$B1

me.es.2 <- all.cd4_es$B3
un.es.2 <- all.cd4_es$B4 - all.cd4_es$B3

dat.mtx = matrix(c(me.cd4.1, un.cd4.1, me.cd4.2, un.cd4.2, me.es.1, un.es.1, me.es.2, un.es.2), ncol = 8)

y = DGEList(dat.mtx)

TotalReadCount <- colMeans(matrix(y$samples$lib.size, nrow=2, ncol=4))
y$samples$lib.size <- rep(TotalReadCount, each=2)

design <- cbind(# Int = 1,
                cd4.2 = c(0,0,1,1,0,0,0,0),
                es.1 = c(0,0,0,0,1,1,0,0),
                es.2 = c(0,0,0,0,0,0,1,1),
                me.cd = c(1,0,1,0,1,0,1,0),
                me.cdVes = c(0,0,0,0,1,0,1,0))

d <- estimateDisp(y, design=design, trend="none")
# y$common.dispersion
d <- estimateGLMCommonDisp(y, design, verbose=TRUE)
fit = glmFit(y, design, dispersion=d)

results <- glmLRT(fit, coef='me.cdVes')
edgeR.pvalue = results$table$PValue

edgeR.qvalue = qvalue(edgeR.pvalue)$qvalues

edgeR.fdr = p.adjust(edgeR.pvalue, method = 'BH')

rej.edgeR = ifelse(edgeR.fdr < 0.1, 1, 0)
rej.edgeR[is.na(rej.edgeR)] = 0 # 981318


#### FDR: 0.05 #####
sum(rej.limma) # 420661
sum(rej.edgeR) # 981318
sum(rej.mcf) # 3959142

sum(rej.limma * rej.mcf) # 412635
sum(rej.edgeR * rej.mcf) # 304382
sum(rej.limma * rej.edgeR) # 77407

#### FDR: 0.1 #####
sum(rej.limma) # 940541
sum(rej.edgeR) # 1073984
sum(rej.mcf) # 5116628

sum(rej.limma * rej.mcf) # 910421
sum(rej.edgeR * rej.mcf) # 406271
sum(rej.limma * rej.edgeR) # 130642





library(ggplot2)

cd4.1 <- all.cd4_es$A2
cd4.2 <- all.cd4_es$A4
cd8.1 <- all.cd8_es$A2
cd8.2 <- all.cd8_es$A4
es.1 <- all.cd4_es$B2
es.2 <- all.cd4_es$B4


dat <- data.frame(sample = factor(c(rep("cd4.1",length(cd4.1)),
                                  rep("cd4.2",length(cd4.2)),
                                  rep("cd8.1",length(cd8.1)),
                                  rep("cd8.2",length(cd8.2)),
                                  rep("es.1",length(es.1)),
                                  rep("es.2",length(es.2)))),
                  total_count = c(cd4.1,cd4.2,cd8.1,cd8.2,es.1,es.2))

p <- ggplot(dat, aes(x=sample, y=total_count)) + 
  geom_boxplot() + 
  coord_cartesian(ylim = c(0, 100))

pdf("/Users/xiaoyudai/Documents/Paper/Real/total_count_box.pdf",width=12,height=8)
p
dev.off()






