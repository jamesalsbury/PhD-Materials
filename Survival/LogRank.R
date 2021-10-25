
library(survival)

survdiff(Surv(time, status) ~ sex, data = lung)


#Need to only include the event times where it is not censored

uniquedead = unique(sort(lung$time[lung$status==2]))
LogRank = data.frame(t = uniquedead, n1 = vector(length=length(uniquedead)), 
           d1 = vector(length=length(uniquedead)),
           n2 = vector(length=length(uniquedead)),
           d2 = vector(length=length(uniquedead)),
           n = vector(length=length(uniquedead)),
           d = vector(length=length(uniquedead)),
           e1 = vector(length=length(uniquedead)),
           d1minuse1 = vector(length=length(uniquedead)),
           V = vector(length=length(uniquedead)))


for (i in 1:length(uniquedead)){
  LogRank[i,]$n1 = sum(lung$time[lung$sex==1]>LogRank[i,]$t)
  LogRank[i,]$d1 = sum(lung$time[lung$status==2&lung$sex==1]==LogRank[i,]$t)
  LogRank[i,]$n2 = sum(lung$time[lung$sex==2]>LogRank[i,]$t)
  LogRank[i,]$d2 = sum(lung$time[lung$status==2&lung$sex==2]==LogRank[i,]$t)
}
LogRank$n = LogRank$n1 + LogRank$n2
LogRank$d = LogRank$d1 + LogRank$d2
LogRank$e1 = (LogRank$n1/LogRank$n)*LogRank$d
LogRank$d1minuse1 = LogRank$d1 - LogRank$e1
LogRank$V = (((LogRank$n1*LogRank$n2)*(LogRank$n-LogRank$d))*LogRank$d)/((LogRank$n^2)*(LogRank$n-1))

T = sum((LogRank$d1minuse1))^2/sum(LogRank$V)
T
