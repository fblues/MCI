sigma<-runif(1000,0.1,10)
hist(sigma)
g<-function(x,b=1)
  x/(b+x)
g(sigma)
gd<-function(x,b=1)
  b/(b+x)^2
ginv<-function(y,b=1)
  b*y/(1-y)

glin<-function(x,h=200,b=1)
{
  n<-length(x)
  cc<-ginv(0:(h-1)/h,b)
  cs<-cut(x,cc)
  labs <- as.vector(cs)
  ds<-cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
            upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
  y<-NULL;slope<-NULL;intercept<-NULL;
  for(i in 1:n){
    slope[i]<-(g(ds[i,2],b)-g(ds[i,1],b))/(ds[i,2]-ds[i,1])
    intercept[i]<--slope[i]*ds[i,1]+g(ds[i,1],b)
    y[i]<-slope[i]*x[i]+intercept[i]
  }
  return(slope)
}

fn_lin<-function(x,sigma,q,b,h)
{
  n<-length(sigma)
  log(2*pnorm(x)-1) + log(glin(sigma*x,h,b)*sigma)+x^2/2 - mean(log(glin(sigma*x,h,b)*sigma))-mean(x^2)/2 - 1/n*log(1-q)
}

mat<-function(x)
{
  n<-length(x)
  diag<-(2*dnorm(x))/(2*pnorm(x)-1)+x
  offdiag<- -(1/n)*x
  mtx<-matrix(rep(offdiag,n),byrow=TRUE,nrow=n)
  for(i in 1:n)
    mtx[i,i]<-diag[i]+mtx[i,i]
  return(mtx)    
}

opt_lin<-function(sigma,b=1,q=0.05,h=200)
{
  n<-length(sigma)
  nu<-rep(qnorm(1-(1-(1-q)^(1/n))/2),n)
  repeat{
    nu_<-matrix(nu)-solve(mat(nu))%*%matrix(fn_lin(nu,sigma,q,b,h))  
    if(sum((nu_-nu)^2)<0.000000001) break
    nu<-nu_}
  eta<-2*(1-pnorm(nu))
  return(cbind(eta,nu))
}


#  beta*se*(2*pnorm(nu)-1)-2*beta*dnorm(nu)*(beta+nu*se)
A<-function(se,nu,beta)
  beta*se

B<-function(se,nu,beta)
  (beta+nu*se)^2
#dA<-function(se,nu,beta)
#  2*beta*nu*dnorm(nu)*(beta+nu*se)
dA<-function(se,nu,beta)
  2*se*dnorm(nu)*(beta+nu*se)*(2-nu^2)

dB<-function(se,nu,beta)
  2*se*(beta+nu*se)



fn_beta<-function(nu,se,beta,q){
  n<-length(nu)
  log(2*pnorm(nu)-1)+log(A(se,nu,beta))-mean(log(A(se,nu,beta)))-log(B(se,nu,beta))+mean(log(B(se,nu,beta)))+nu^2/2-mean(nu^2/2)-1/n*log(1-q)
}

mat_beta<-function(nu,se,beta){
  n<-length(nu)
  diag<-(2*dnorm(nu))/(2*pnorm(nu)-1)-dB(se,nu,beta)/B(se,nu,beta)+nu
  offdiag<--(1/n)*(-dB(se,nu,beta)/B(se,nu,beta)+nu)
  mtx<-matrix(rep(offdiag,n),byrow=TRUE,nrow=n)
  for(i in 1:n)
    mtx[i,i]<-diag[i]+mtx[i,i]
  return(mtx)    
}

opt_beta<-function(se,beta=1,q=0.05){
  n<-length(se)
  nu<-rep(qnorm(1-(1-(1-q)^(1/n))/2),n)
  repeat{
    nu_new<-matrix(nu)-solve(mat_beta(nu,se,beta))%*%matrix(fn_beta(nu,se,beta,q))  
    if(sum((nu_new-nu)^2)<0.0000001) break
    nu<-nu_new}
  eta<-2*(1-pnorm(nu))
  return(cbind(eta,nu))
}


pic_b<-function(b=c(1,5,10),sigma=sigma){
  n<-length(b); m<-length(sigma)
  eta<-matrix(rep(0,m*(n+1)),ncol=n+1);nu<-matrix(rep(0,m*(n+1)),ncol=n+1)
  eta_lin<-matrix(rep(0,m*(n+1)),ncol=n+1);nu_lin<-matrix(rep(0,m*(n+1)),ncol=n+1)
  for(i in 1:n)
  {
    out_lin<-opt_lin(sigma,b[i],h=500)
    out<-opt_beta(sigma,b[i])
    eta_lin[,i]<-out_lin[,1];nu_lin[,i]<-out_lin[,2]
    eta[,i]<-out[,1];nu[,i]<-out[,2]
  }
  nu[,n+1]<-rep(qnorm(1-(1-(1-0.05)^(1/m))/2),m)
  eta[,n+1]<-rep((1-(1-0.05)^(1/m)),m)
  nu_lin[,n+1]<-rep(qnorm(1-(1-(1-0.05)^(1/m))/2),m)
  eta_lin[,n+1]<-rep((1-(1-0.05)^(1/m)),m)
  return(list(eta_lin,eta,nu_lin,nu))
}
pic_h<-function(b=1000,sigma=sigma,h=c(50,100,500)){
  n<-length(h); m<-length(sigma)
  eta<-matrix(rep(0,m*(n+1)),ncol=n+1);nu<-matrix(rep(0,m*(n+1)),ncol=n+1)
  eta_lin<-matrix(rep(0,m*(n+1)),ncol=n+1);nu_lin<-matrix(rep(0,m*(n+1)),ncol=n+1)
  for(i in 1:n)
  {
    out_lin<-opt_lin(b=b,sigma,h[i])
    eta_lin[,i]<-out_lin[,1];nu_lin[,i]<-out_lin[,2]
  }
  eta_lin[,n+1]<-opt_beta(sigma,b)[,1]
  nu_lin[,n+1]<-opt_beta(sigma,b)[,2]
  return(list(eta_lin,nu_lin,out))
}
opt_h50<-opt_lin(b=2,sigma,h=50)
opt_h100<-opt_lin(b=2,sigma,h=100)
opt_h500<-opt_lin(b=2,sigma,h=500)
opt_h1000<-opt_beta(b=2,sigma,q=0.05)
matplot(sigma,data.frame(opt_h50[,1],opt_h100[,1],opt_h500[,1],opt_h1000[,1],opt1_h50[,1],opt1_h100[,1],opt1_h500[,1],opt1_h1000[,1]),type="p",col=c("purple","blue","red","black","purple","blue","red","black"),pch=19,cex=0.5,ann=FALSE)
legend("topright",legend=c(expression(gamma==50),expression(gamma==100),expression(gamma==500),"Original Ftn"),col=c("purple","blue","red","black"),lwd=3)
title(main="Interpolation Quality",ylab="Size",xlab="Sigma",cex.lab=1.25,cex.main=1.7)

aa<-pic(c(2,5,10,50,1000),sigma)
bb<-pic_h(b=1000,sigma,c(100,200,500))
matplot(sigma,data.frame(aa[1]),type="p",col=c("black","green","blue","purple","red","grey"),pch=19,cex=0.7,ann=FALSE)
legend("topright",legend=c(expression(beta==2),expression(beta==5),expression(beta==10),expression(beta==50),expression(beta==1000),"Sidak"),,col=c("black","green","blue","purple","red","grey"),lwd=3)
title(main="Size Investing Strategy",ylab="Size",xlab="Sigma",cex.lab=1.25,cex.main=1.8)

matplot(sigma,data.frame(aa[2]),type="p",col=c("black","green","blue","purple","red","grey"),pch=19,cex=0.7,main="Size Allocation: Original Loss Function",ylab="Size")
legend("topright",legend=c(expression(beta==2),expression(beta==6),expression(beta==10),expression(beta==30),expression(beta==300),"Sidak"),,col=c("black","green","blue","purple","red","grey"),lwd=3)
matplot(sigma,data.frame(aa[3]),type="p",col=c("black","green","blue","purple","red","grey"),pch=19,cex=0.7,ann=FALSE)
legend("topright",legend=c(expression(beta==2),expression(beta==6),expression(beta==10),expression(beta==30),expression(beta==300),"Sidak"),,col=c("black","green","blue","purple","red","grey"),lwd=3)
title(main="Confidence Coeffefficient",ylab="Confidence Coeff.",xlab="Sigma",cex.lab=1.25,cex.main=1.7)


matplot(sigma,data.frame(aa[4]),type="p",col=c("black","green","blue","purple","red","grey"),pch=19,cex=0.7,main="Confidence Coeff.: Original Loss Function",ylab="Confidence Coeff.")
legend("topright",legend=c(expression(beta==2),expression(beta==6),expression(beta==10),expression(beta==30),expression(beta==300),"Sidak"),,col=c("black","green","blue","purple","red","grey"),lwd=3)

matplot(sigma,data.frame(aa[3])*sigma,type="p",ylim=c(20,45),xlim=c(5,10),col=c("black","green","blue","purple","red","grey"),pch=19,cex=0.7,ann=FALSE)
legend("bottomright",legend=c(expression(beta==2),expression(beta==5),expression(beta==10),expression(beta==50),expression(beta==1000),"Sidak"),,col=c("black","green","blue","purple","red","grey"),lwd=3)
title(main="Individual Interval Length",ylab="Length",xlab="Sigma",cex.lab=1.25,cex.main=1.8)
colSums(data.frame(aa[3])*sigma)
M=100
M
para
k<-Out(100,b=c(1,2),q=0.05,para,h=c(50,100))
Out<-function(M,b,q,par,h)
{
  # Data Generation
  Data<-Gen_Data_normal(M,par)
  # NR optimization
  n<-length(b)
  m<-length(h)
#  Output<-NULL;#list(mode="list",m+1);CI<-vector(mode="list",m+1)
 # Length<-NULL;#list();Inclusion<-list();#In<-NULL;Len<-NULL
  #In<-NULL;Len<-NULL
  nu_Sidak<-qnorm(1-(1-(1-q)^(1/M))/2)
  In<-matrix(rep(0,(n+1)*(m+1)),ncol=n+1)
  Len<-matrix(rep(0,(n+1)*(m+1)),ncol=n+1)
  for(j in 1:m){
  for(i in 1:n){
    Output<-opt_lin(sigma=Data[,2],q=q,b=b[i],h=h[j])
    CI<-cbind(Data[,3]-Output[,2]*Data[,2],Data[,3]+Output[,2]*Data[,2])
    Length<-cbind(CI[,2]-CI[,1])
    Inclusion<-cbind(CI[,1]<=Data[,1] & CI[,2]>=Data[,1])
    In[j,i]<-(sum(Inclusion)==M)
    Len[j,i]<-mean(Length)
  }
  }
  for(i in 1:n){
  Output<-opt_beta(Data[,2],b[i],q=q)
  CI<-cbind(Data[,3]-Output[,2]*Data[,2],Data[,3]+Output[,2]*Data[,2])
  Length<-cbind(CI[,2]-CI[,1])
  Inclusion<-cbind(CI[,1]<=Data[,1] & CI[,2]>=Data[,1])
  In[m+1,i]<-(sum(Inclusion)==M)
  Len[m+1,i]<-mean(Length)
  }
  In;Len
  CI<-cbind(Data[,3]-nu_Sidak*Data[,2],Data[,3]+nu_Sidak*Data[,2])
  Length<-CI[,2]-CI[,1]
  Inclusion<-cbind(CI[,1]<=Data[,1] & CI[,2]>=Data[,1])
  In[m+1,n+1]<-(sum(Inclusion)==M)
  Len[m+1,n+1]<-mean(Length)
  return(list(In,Len))
}

Gen_Data_normal<-function(M,par)
{
  dat_list<-NULL  
  for(i in 1:M)
    dat_list[i]<-rnorm(1,par[i,1],par[i,2])
  return(cbind(par[,1],par[,2],dat_list))  
}

1-(1-0.05)^(1/2000)

Final<-function(re=1000,M=51,b=c(2,50,1000),q=0.05,h=c(50,100,500))
{
  mu_vec<-rnorm(M,0,3);sig_vec<-runif(M,0.01,10)
  para<-cbind(mu_vec,sig_vec)
  n<-(length(b)+1)
  m<-(length(h)+1)
  In<-matrix(rep(0,n*m*re),ncol=n*m)
  Len<-matrix(rep(0,n*m*re),ncol=n*m)
  for(i in 1:re){
    cat("loop #",i,"\n")
    output<-Out(M,b,q,par=para,h)
    In[i,]<-as.vector(output[[1]])
    Len[i,]<-as.vector(output[[2]])
  }
  mn_In<-apply(In,2,mean)
  mn_Len<-apply(Len,2,mean)
  sd_In<-apply(In,2,sd)
  sd_Len<-apply(Len,2,sd)
  mnIN<-matrix(mn_In,ncol=n);mnLN<-matrix(mn_Len,ncol=n)
  sdIN<-matrix(sd_In,ncol=n);sdLN<-matrix(sd_Len,ncol=n)
   # Label<-c("Cover_beta1:","TLength_beta1:","Cover_beta2:","TLength_beta2:","Cover_beta3:","TLength_beta3:","Cover_test:","TLength_test:","Cover_Sidak:","TLength_Sidak:")
  return(list(mnIN,mnLN,sdIN,sdLN))
}
Final(re=10,b=c(2,10,1000),h=500,M=100)

Final(re=1000,b=c(2,10,1000),q=0.05,h=500,M=500)
matrix(as.vector(k[[1]]),ncol=3)
