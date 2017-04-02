load("C:/Dropbox/Public/Research/Research/Bootstrap/rate_of_return.Rdata")

return

a6<-boot_CI(return,,B=5000,b=50)

boot_CI<-function(data,q=0.05,B=1000,b=20)
{
  id<-dim(data)
  CIs<-matrix(rep(0,4*id[2]),ncol=4)
  for(i in 1:id[2])
  {
  print(i)
  CIs[i,]<-Ind_CI(data[,i],q,B,b)
  }
  return(CIs)
}




Ind_CI<-function(company,q,B,b)
{
med<-median(company)
n<-length(company)
medB<-NULL;zboot<-NULL;medBb<-matrix(rep(0,B*b),ncol=b)
for(j in 1:B)
{
  sample_B<-sample(company,n,replace=TRUE)
  medB[j]<-median(sample_B)
  for(k in 1:b)
  {
    sample_b<-sample(sample_B,n,replace=TRUE)
    medBb[j,k]<-median(sample_b)
  }
  zboot[j]<-(medB[j]-med)/sd(medBb[j,])
}
CC<-as.numeric(quantile(zboot,c(q/2,1-q/2),na.rm=TRUE))
SE<-sd(medB)
CI<-c(med,CC,SE)
return(CI)
}


cbind(a1[,4],a2[,4])
sum(abs(a1[,2]-a2[,2]))
sum(abs(a3[,2]-a4[,2]))
sum(abs(a1[,4]-a2[,4]))
sum(abs(a3[,4]-a4[,4]))
sum(abs(a5[,4]-a6[,4]))
a6
a5
cbind(a5[,4],a6[,4])
