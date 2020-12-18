library(survival)
library(ggplot2)
library(survminer)
library(grid)
library(mvtnorm)
library(stats4)
library(modeltools)
library(zoo)
library(strucchange)
library(sandwich)
library(party)
library(maxstat)

#generate OS of COPS5
setwd("~/jab1_cops5")
methy<-read.table("OS_COPS5.txt",header=TRUE)
stat <- maxstat.test(Surv(OS, OS_IND)~cg10262052, data=methy, smethod="LogRank", pmethod="exactGauss", abseps=0.01)
stat
plot(stat)#出来cutpoint， 代入下面一行的值
methygroup = ifelse(methy$cg10262052<=10.237, 0, 1)
methygroup = factor(methygroup)
levels(methygroup) = c("below 10.237", "above 10.237")
fit <- survfit(Surv(OS, OS_IND) ~ methygroup, data=methy)
fit
plot(fit)
survdiff(Surv(OS, OS_IND)~methygroup, data=methy)
out=coxph(Surv(OS, OS_IND)~methygroup, data=methy)
out
ggsurvplot(fit, pval=TRUE, legent.title="cg10262052",
           legend.labs=c("below cut-point", "above cut-point"),
           xlab="Time(days)", break.time.by=1000, main="TCGA LAML hMethyl 450 - cg10262052",
           legend="right", palette=c("#00BFC4", "#000000"))


#generate RFS of COPS5
methy<-read.table("RFS_COPS5.txt",header=TRUE)
stat <- maxstat.test(Surv(OS, OS_IND)~cg10262052, data=methy, smethod="LogRank", pmethod="exactGauss", abseps=0.01)
stat
plot(stat)#出来cutpoint， 代入下面一行的值
methygroup = ifelse(methy$cg10262052<=10.796, 0, 1)
methygroup = factor(methygroup)
levels(methygroup) = c("below 10.796", "above 10.796")
fit <- survfit(Surv(OS, OS_IND) ~ methygroup, data=methy)
fit
plot(fit)
survdiff(Surv(OS, OS_IND)~methygroup, data=methy)
out=coxph(Surv(OS, OS_IND)~methygroup, data=methy)
out
ggsurvplot(fit, pval=TRUE, legent.title="cg10262052",
           legend.labs=c("below cut-point", "above cut-point"),
           xlab="Time(days)", break.time.by=1000, main="TCGA LAML hMethyl 450 - cg10262052",
           legend="right", palette=c("#00BFC4", "#000000"))