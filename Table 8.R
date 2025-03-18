require(randomForest)
require(ROCR)

require(cluster)
require(rlist)
require(ggplot2)
require(plotly)
require(philentropy)
library(MASS)
library(clusterGeneration)
library(clusterCrit)
library(clusterSim)
library(clValid)
library(GIC)
library(MVTests)
library(vegan)
library(KmeansInference)
library(readr)
library(mclust)
library(VALIDICLUST)
library(ggpubr)
library(DescTools)
library(CADET)

clus = function(km,newdata){
  max = km$iter
  c=km$centers[[max]]
  
  aa=matrix(NA,nrow(newdata),nrow(c))
  
  for(i in 1:nrow(c)){
    uk=data.frame(t(c[i,]))
    uk=uk[rep(seq_len(nrow(uk)), each = nrow(newdata)),]
    aa[,i]=sqrt(rowSums((newdata-uk)^2))
  }
  bb=aa==apply(aa,1,min)
  cc=NA
  for(i in 1:nrow(newdata)){
    cc[i]=c(1:nrow(c))[bb[i,]]
  }
  
  result=list("dis"=km$final_cluster,"val"=cc)
  result
}

##Flip
flip =function(x){
  x=x+1
  x=ifelse(x==3,1,x)
  x
}



sg=matrix(c(1,-0.7,0.3,0.131,
            -0.7,1,0.05,0.59,
            0.3,0.05,1,0.34,
            0.131,0.59,0.34,1),
          nrow=4,ncol=4)


mu2= rep(3,4)
mu3= rep(-3,4)

rra=matrix(0,9,4)
rrb=matrix(0,9,5)

rsa = matrix(0,1,3)
rsb = matrix(0,1,4)

for(i in 1:1000){
  
  
  x2=mvrnorm(n = 100, mu=mu2, Sigma =sg)
  x3=mvrnorm(n = 100, mu=mu3, Sigma =sg)
  d=as.matrix(rbind(x2,x3))
  
  res_Kmeans = KmeansInference::kmeans_estimation(as.matrix(d), 2, seed = NULL)
  cl1 = res_Kmeans$final_cluster
  za=aggregate(d[,1],list(cl1),mean)[1,2]
  if(za>0){
    cl1 = flip(cl1)
  }
  
  
  ##D1
  
  
  uuu=kmeans_inference(d[cl1==1,],2,1,2,iso=F,SigInv =solve(cov(d[cl1==1,])),seed=NULL)$pval
  
  rsa[1,1]=rsa[1,1]+sum(uuu>0.05)
  
  ##D2
  
  
  uuu=kmeans_inference(d[cl1==2,],2,1,2,iso=F,SigInv =solve(cov(d[cl1==2,])),seed=NULL)$pval
  
  rsa[1,2]=rsa[1,2]+sum(uuu>0.05)
  
  ## D1 vs D2
  
  
  
  uuu  =kmeans_inference(d,2,1,2,iso=F,SigInv =solve(cov(d)),seed=NULL)$pval
  
  rsa[1,3]=rsa[1,3]+sum(uuu>0.05)
  
  
  for(ii in 1:4){
    
    
    g1=kmeans_inference_1f(as.matrix(d),k = 2, 1,2,
                           feat = ii, iso = F,
                           covMat=cov(d),
                           iter.max = 15,seed=NULL)$pval
    
    
    
    rsb[1,ii]=rsb[1,ii]+sum(g1>0.05)
    
  }
  
  
  uq=1
  for(k in seq(1,9,1)){
    rra[uq,1]=k
    rrb[uq,1]=k
    
    mu22=rep(3,4)
    mu33= rep(-k,4)
    
    ######Validation
    
    x2a=mvrnorm(n = 100, mu=mu22, Sigma =sg)
    x3a=mvrnorm(n = 100, mu=mu33, Sigma =sg)
    v=as.matrix(rbind(x2a,x3a))
    
    zz=clus(res_Kmeans,v)
    cl2=zz$val
    
    if(za>0){
      cl2 = flip(cl2)
    }
    
    xx=rbind(d,v)
    ct = c(cl1,cl2)
    
    
    uuu   =kmeans_inference(xx[ct==1,],2,1,2,iso=F,SigInv =solve(cov(xx[ct==1,])),seed=NULL)$pval
    
    rra[uq,2]=rra[uq,2]+sum(uuu>0.05)
    
    ##D2 U V2
    
    
    
    uuu   =kmeans_inference(xx[ct==2,],2,1,2,iso=F,SigInv =solve(cov(xx[ct==2,])),seed=NULL)$pval
    
    rra[uq,3]=rra[uq,3]+sum(uuu>0.05)
    
    ##D1 V1 U D2 V2
    
    
    uuu  = kmeans_inference(xx,2,1,2,iso=F,SigInv =solve(cov(xx)),seed=NULL)$pval
    
    rra[uq,4]=rra[uq,4]+sum(uuu>0.05)
    
    
    for(ii in 1:4){
      
      
      
      g2=kmeans_inference_1f(as.matrix(xx),k = 2, 1,2,
                             feat = ii, iso = F,
                             covMat=cov(xx),
                             iter.max = 15,seed=NULL)$pval
      
      
      
      
      rrb[uq,ii+1]=rrb[uq,ii+1]+sum(g2>0.05)
      
      
    }
    
    
    
    uq=uq+1
  }
  
}
write.csv(rrb/1000,"D:\\Working\\NYU Project\\Carole Validation\\New test 20240904\\New Paper Tables\\final_0207\\gamma_feature.csv",row.names = F)
