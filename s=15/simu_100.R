set.seed( 100 ) 
  library(mvtnorm)
  library(fields)
  library(Rcpp)
  library(mclust)
  library(kernlab)
  library(ConsensusClusterPlus)
  
  simu=function(s){
  
  prob_glcm<-function(c,s=s,mc=30000){
  mu<-c(2+c,14-c)
  sigma<-matrix(s*c(1,-0.7,-0.7,1),nrow=2)
  elip<-rmvnorm(mc,mu,sigma)
  # par(xaxs='i',yaxs='i')
  # plot(elip,xlim =c(0,16) ,ylim=c(0,16))
  # abline(16,-1,col='red')
  # abline(h=16);abline(h=15);abline(h=14);abline(h=13);abline(h=12);abline(h=11);abline(h=10);abline(h=9);
  # abline(h=8);abline(h=7);abline(h=6);abline(h=5);abline(h=4);abline(h=3);abline(h=2);abline(h=1);abline(h=0)
  # abline(v=16);abline(v=15);abline(v=14);abline(v=13);abline(v=12);abline(v=11);abline(v=10);abline(v=9);
  # abline(v=0);abline(v=1);abline(v=2);abline(v=3);abline(v=4);abline(v=5);abline(v=6);abline(v=7);abline(v=8)
  
  cell_count<-rep(0,16*16)
  
  for (i in 1:mc)
  {
  for (m in 1:16) {
  
  for (k in 16:1) {
  
  if (( (m-1) <elip[i,1])&(elip[i,1]< m)&( (k-1) <elip[i,2])&(elip[i,2]< k)) {
  cell_count[16-k+1+16*(m-1)]=cell_count[16-k+1+16*(m-1)]+1}
  }
  
  }
  }
  
  ## -c(2:16,19:32,36:48,53:64,70:80,87:96,104:112,121:128,138:144,155:160,172:176,189:192,206:208,223:224,240)
  
  
  z<-cell_count/sum(cell_count)
  z_whole<-z[c(1,17,33,49,65,81,97,113,129,145,161,177,193,209,225,241,
  17,18,34,50,66,82,98,114,130,146,162,178,194,210,226,242,
  33,34,35,51,67,83,99,115,131,147,163,179,195,211,227,243,
  49,50,51,52,68,84,100,116,132,148,164,180,196,212,228,244,
  65,66,67,68,69,85,101,117,133,149,165,181,197,213,229,245,
  81,82,83,84,85,86,102,118,134,150,166,182,198,214,230,246,
  97,98,99,100,101,102,103,119,135,151,167,183,199,215,231,247,
  113,114,115,116,117,118,119,120,136,152,168,184,200,216,232,248,
  129,130,131,132,133,134,135,136,137,153,169,185,201,217,233,249,
  145,146,147,148,149,150,151,152,153,154,170,186,202,218,234,250,
  161,162,163,164,165,166,167,168,169,170,171,187,203,219,235,251,
  177,178,179,180,181,182,183,184,185,186,187,188,204,220,236,252,
  193,194,195,196,197,198,199,200,201,202,203,204,205,221,237,253,
  209,210,211,212,213,214,215,216,217,218,219,220,221,222,238,254,
  225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,255,
  241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256)]
  arg. <- expand.grid(c(0.5:15.5),c(15.5:0.5))
  I = as.image( Z=z_whole, x=arg., grid=list(x=seq(0.5,15.5,1), y=seq(0.5,15.5,1)))
  image(I)
  smooth.I <- image.smooth(I, theta=1);
  
  
  #################################################
  ###  notice the order of this sommthed image  ###
  #################################################  
  den=c()
  for (r in 1:16) {
  for (w in 1:r) {
  den=c(den,smooth.I$z[r,16-(w-1)])
  }
  }
  
  
  prob<-den/sum(den)
  return(prob) 
  }
  
  
  prob1=prob_glcm(c=5,s=s)
  prob2=prob_glcm(c=5.5,s=s)
  prob3=prob_glcm(c=6,s=s)
  prob4=prob_glcm(c=6.5,s=s)
  prob5=prob_glcm(c=7,s=s)

  glcm=matrix(0,nrow=20*5,ncol=136)
  for (j in 1:20)
  {
  t<-round(runif(1,500,20000),0)
  glcm[j,]=round(t*prob1)
  }
  
  for (j in 21:40)
  {
  t<-round(runif(1,500,20000),0)
  glcm[j,]=round(t*prob2)
  }
  
  for (j in 41:60)
  {
  t<-round(runif(1,500,20000),0)
  glcm[j,]=round(t*prob3)
  }
  
  for (j in 61:80)
  {
  t<-round(runif(1,500,20000),0)
  glcm[j,]=round(t*prob4)
  } 
  
  for (j in 81:100)
  {
  t<-round(runif(1,500,20000),0)
  glcm[j,]=round(t*prob5)
  }
  
  glcm
  
  }
  
  Z=simu(s=15)

  Z_met=Z
  T_met=nrow(Z_met)
  n=ncol(Z_met)
  
  X=apply(Z_met,1,sum)
  X_met=X
  sX_met=(X-mean(X))/sd(X)
  
  R=array(data = NA,dim = c(2,n,T_met))
  for (t in 1: nrow(Z_met)) R[,,t]=matrix(rep(c(1,sX_met[t]),times=n),byrow = FALSE,nrow=2,ncol=n)
  
  
  ############################################################################
  ##########################            MCMC          ########################
  ############################################################################
  
  
  library(HI)
  library(invgamma)
  
  source('/gstore/scratch/u/lix233/RGSDP/sdp_functions_selfwriting_V12_cpp.R')
  sourceCpp('/gstore/scratch/u/lix233/RGSDP/rgsdp.cpp')
  
  D=read.csv('/gstore/scratch/u/lix233/RGSDP/D_16.csv',header=TRUE)
  W=read.csv('/gstore/scratch/u/lix233/RGSDP/W_16.csv',header=TRUE)
  
  N=20000;Burnin=N/2
  
  Y_iter_met=Theta_iter_met=array(data=NA,dim = c(T_met,n,N))
  
  
  try=matrix(0,nrow =T_met ,ncol = n)
  
  for (i in 1:T_met){
  for (j in 1:n){
  if (Z_met[i,j]==0) {
  try[i,j]=rnorm(1,mean=-10,sd=1)
  } else {
  try[i,j]=rnorm(1,mean=Z_met[i,j],sd=1)
  }
  }
  } 
  
  
  
  g=update_Y(Z=Z_met,X=X_met,tau2=100,Theta = try,Beta =c(0.1,0.1),R)
  sum(g==Inf)+sum(g==-Inf)
  
  
  Theta_iter_met[,,1]=try
  
  tau2_met=v_met=rho_met=sig2_met=rep(NA,N)
  tau2_met[1]=50
  v_met[1]=0.8
  rho_met[1]=0.9
  sig2_met[1]=10
  
  # v_met=rep(1,N) # Fix v
  

  av=bv=1
  atau=0.0001 ;btau=0.0001
  asig=0.0001 ;bsig=0.0001
  
  Betam=c(0,0);Sigma_m=matrix(c(10^5,0,0,10^5),nrow=2,ncol=2)
  
  Beta_iter_met=matrix(NA,nrow=N,ncol=nrow(R[,,1]))
  
  Beta_iter_met[1,]=c(40,20)
  
  
  
  
  
  
  
  for (iter in 2:N) {
  
  
  Y_iter_met[,,iter]=update_Y(Z_met,X_met,tau2_met[iter-1],Theta_iter_met[,,iter-1],Beta_iter_met[iter-1,],R)
  
  Theta_iter_met[,,iter]=update_theta(as.vector(X_met),Y_iter_met[,,iter],as.matrix(D),as.matrix(W),rho_met[iter-1],Theta_iter_met[,,iter-1],sig2_met[iter-1],tau2_met[iter-1],v_met[iter-1],Beta_iter_met[iter-1,],R)
  
  Beta_iter_met[iter,]=update_Beta(Betam,Sigma_m,tau2_met[iter-1],X_met,Y_iter_met[,,iter],Theta_iter_met[,,iter],R)
  
  tau2_met[iter] = update_tau2(X_met,Y_iter_met[,,iter],Theta_iter_met[,,iter],atau,btau,Beta_iter_met[iter,],R)
  
  sig2_met[iter]= update_sig2(asig,bsig,D,W,rho_met[iter-1],Theta_iter_met[,,iter])
  
  rho_met[iter] = update_rho(D,W,Theta_iter_met[,,iter],sig2_met[iter])
  
  v_met[iter]=update_v(Z_met,v_met[iter-1],Tstar=nrow(unique.matrix(Theta_iter_met[,,iter])),av,bv)

 
  
  }
  
  
  library(coda)
  
  mcmc_beta=mcmc(Beta_iter_met[(1+Burnin):N,])
  pnorm(abs(geweke.diag(mcmc_beta)$z),lower.tail=FALSE)*2
  
  
  
  mcmc_rho=mcmc(rho_met[(1+Burnin):N])
  pnorm(abs(geweke.diag(mcmc_rho)$z),lower.tail=FALSE)*2
  
  
  
  mcmc_sig2=mcmc(sig2_met[(1+Burnin):N])
  pnorm(abs(geweke.diag(mcmc_sig2)$z),lower.tail=FALSE)*2
  
  
  mcmc_tau2=mcmc(tau2_met[(1+Burnin):N])
  pnorm(abs(geweke.diag(mcmc_tau2)$z),lower.tail=FALSE)*2
 

  mcmc_v=mcmc(v_met[(1+Burnin):N])
  pnorm(abs(geweke.diag(mcmc_v)$z),lower.tail=FALSE)*2  

  
  
  
  Theta_ave=Theta_sum=matrix(0,nrow=nrow(Theta_iter_met[,,1]),ncol=ncol(Theta_iter_met[,,1]))
  for (i in (Burnin+1):N) {
  Theta_sum=Theta_sum+Theta_iter_met[,,i]
  }
  
  Theta_ave=Theta_sum/(N-Burnin)
  
  
 
  
  library('NbClust')
  NbClust(Theta_ave,distance='euclidean',method='ward.D2',index='kl')
  
  HRGSDP=NbClust(Theta_ave,distance='euclidean',method='ward.D2',index='kl')$Best.partition
  
  
  
  
  
  glcm_whole=Z[,c(1,2,4,7,11,16,22,29,37,46,56,67,79,92,106,121,
  2,3,5,8,12,17,23,30,38,47,57,68,80,93,107,122,
  4,5,6,9,13,18,24,31,39,48,58,69,81,94,108,123,
  7,8,9,10,14,19,25,32,40,49,59,70,82,95,109,124,
  11,12,13,14,15,20,26,33,41,50,60,71,83,96,110,125,
  16,17,18,19,20,21,27,34,42,51,61,72,84,97,111,126,
  22,23,24,25,26,27,28,35,43,52,62,73,85,98,112,127,
  29,30,31,32,33,34,35,36,44,53,63,74,86,99,113,128,
  37,38,39,40,41,42,43,44,45,54,64,75,87,100,114,129,
  46,47,48,49,50,51,52,53,54,55,65,76,88,101,115,130,
  56,57,58,59,60,61,62,63,64,65,66,77,89,102,116,131,
  67,68,69,70,71,72,73,74,75,76,77,78,90,103,117,132,
  79,80,81,82,83,84,85,86,87,88,89,90,91,104,118,133,
  92,93,94,95,96,97,98,99,100,101,102,103,104,105,119,134,
  106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,135,
  121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136)]
  
  
  source('/gstore/scratch/u/lix233/RGSDP/cal_stat.R')
  
  features=cal_stat(glcm_whole)
  
  GMM=Mclust(features,5)
  my.dist <- function(x) dist(x, method='euclidean')
  my.hclust <- function(d) hclust(d, method='ward.D2')
  HC<-cutree(my.hclust(my.dist(data.matrix(features))),k=5)
  KM=kmeans(features,5)
  SC=specc(features,5)
  CO <- ConsensusClusterPlus(t(features),maxK=9,reps=100,pItem=0.90, pFeature=1,
  clusterAlg='hc',distance='euclidean',plot=FALSE)
  CO <- CO[[5]]$consensusClass
  
  
  aa <- table(rep(1:5,each=20), CO)
  bb <- table(rep(1:5,each=20), GMM$classification)  
  cc <- table(rep(1:5,each=20), HC)
  dd <- table(rep(1:5,each=20), KM$cluster)
  ee <- table(rep(1:5,each=20), SC)  
  ff <- table(rep(1:5,each=20), HRGSDP)
  
  
  res_FeaCO=c(chisq.test(aa,correct = TRUE)$statistic,ncol(aa),error_rate(aa), 'FeaCO')
  res_FeaGMM=c(chisq.test(bb,correct = TRUE)$statistic,ncol(bb),error_rate(bb), 'FeaGMM')
  res_FeaHC=c(chisq.test(cc,correct = TRUE)$statistic,ncol(cc),error_rate(cc), 'FeaHC')
  res_FeaKM=c(chisq.test(dd,correct = TRUE)$statistic,ncol(dd),error_rate(dd), 'FeaKM')
  res_FeaSC=c(chisq.test(ee,correct = TRUE)$statistic,ncol(ee),error_rate(ee), 'FeaSC')
  res_HRGSDP=c(chisq.test(ff,correct = TRUE)$statistic,ncol(ff),error_rate(ff), 'HRGSDP')
  
  
  
  xx = rbind(res_FeaCO, res_FeaGMM, res_FeaHC, res_FeaKM, res_FeaSC, res_HRGSDP)
  
  colnames(xx) = c('pearson.chi.sq', 'nunber of clusters', 'error.rate', 'method')

  xx = as.data.frame(xx)
  print(xx)
