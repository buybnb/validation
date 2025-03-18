

####For Three Clusters

library(KmeansInference)
library(palmerpenguins)
library(fpc)
library(CADET)
t=penguins
t=t[,c(1,3:6)]
t=t[complete.cases(t),]
ta=t$species
t=t[,-1]
t=scale(t)
set.seed(1985)
ii=sample(1:nrow(t),200)
t1=t[ii,]
t2=t[-ii,]
ua1=ta[ii]
ub1=ta[-ii]

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

###

tu=kmeans_estimation(as.matrix(t1), 3,iter.max = 10,seed=1234)
table(tu$final_cluster,ua1)

uz=clus(tu,t2)
table(uz$val,ub1)

##D1
#ua=tu$centers[[3]][c(2,1),]
#sqrt(sum((ua[1,]-ua[2,])^2))
t1=data.frame(t1)
t1$ct=drop(tu$final_cluster)

te1=t1[t1$ct==1,]
aaa=kmeans(te1[,1:4],2)
###Distance
sqrt(sum(((aaa$centers[1,])-(aaa$centers[2,]))^2))

##p-value

kmeans_inference(as.matrix(te1[,1:4]),2,1,2,iso=F,SigInv =solve(cov(as.matrix(te1[,1:4])) ))$pval

##D2
#ua=tu$centers[[3]][c(2,1),]
#sqrt(sum((ua[1,]-ua[2,])^2))


te1=t1[t1$ct==2,]
aaa=kmeans(te1[,1:4],2)
###Distance
sqrt(sum(((aaa$centers[1,])-(aaa$centers[2,]))^2))

##p-value

kmeans_inference(as.matrix(te1[,1:4]),2,1,2,iso=F,SigInv =solve(cov(as.matrix(te1[,1:4])) ))$pval


##D3
#ua=tu$centers[[3]][c(2,1),]
#sqrt(sum((ua[1,]-ua[2,])^2))


te1=t1[t1$ct==3,]
aaa=kmeans(te1[,1:4],2)
###Distance
sqrt(sum(((aaa$centers[1,])-(aaa$centers[2,]))^2))


kmeans_inference(as.matrix(te1[,1:4]),2,1,2,iso=F,SigInv =solve(cov(as.matrix(te1[,1:4])) ))$pval

##
ua=tu$centers[[tu$iter]][c(2,1),]
sqrt(sum((ua[1,]-ua[2,])^2))


kmeans_inference(as.matrix(t1[,1:4]),3,1,2,iso=F,SigInv =solve(cov(as.matrix(t1[,1:4])) ))$pval

  
  
ua=tu$centers[[tu$iter]][c(3,1),]
sqrt(sum((ua[1,]-ua[2,])^2))

kmeans_inference(as.matrix(t1[,1:4]),3,1,3,iso=F,SigInv =solve(cov(as.matrix(t1[,1:4])) ))$pval


ua=tu$centers[[tu$iter]][c(3,2),]
sqrt(sum((ua[1,]-ua[2,])^2))

kmeans_inference(as.matrix(t1[,1:4]),3,3,2,iso=F,SigInv =solve(cov(as.matrix(t1[,1:4])) ))$pval

##
t2=data.frame(t2)
t2$ct=uz$val

tt = rbind(t1,t2)


te1=tt[tt$ct==1,]
aaa=kmeans(te1[,1:4],2)
###Distance
sqrt(sum(((aaa$centers[1,])-(aaa$centers[2,]))^2))

##p-value
kmeans_inference(as.matrix(te1[,1:4]),2,1,2,iso=F,SigInv =cov(as.matrix(te1[,1:4])) )

te1=tt[tt$ct==2,]
aaa=kmeans(te1[,1:4],2)
###Distance
sqrt(sum(((aaa$centers[1,])-(aaa$centers[2,]))^2))

##p-value
kmeans_inference(as.matrix(te1[,1:4]),2,1,2,iso=F,SigInv =cov(as.matrix(te1[,1:4])) )


te1=tt[tt$ct==3,]
aaa=kmeans(te1[,1:4],2)
###Distance
sqrt(sum(((aaa$centers[1,])-(aaa$centers[2,]))^2))

##p-value
kmeans_inference(as.matrix(te1[,1:4]),2,1,2,iso=F,SigInv =cov(as.matrix(te1[,1:4])) )


####
sqrt(sum((apply(tt[tt$ct==1,1:4],2,mean)-apply(tt[tt$ct==2,1:4],2,mean))^2))

kmeans_inference(as.matrix(tt[,1:4]),3,1,2,iso=F,SigInv =cov(as.matrix(tt[,1:4])) )


sqrt(sum((apply(tt[tt$ct==1,1:4],2,mean)-apply(tt[tt$ct==3,1:4],2,mean))^2))

kmeans_inference(as.matrix(tt[,1:4]),3,1,3,iso=F,SigInv =cov(as.matrix(tt[,1:4])) )

sqrt(sum((apply(tt[tt$ct==2,1:4],2,mean)-apply(tt[tt$ct==3,1:4],2,mean))^2))

kmeans_inference(as.matrix(tt[,1:4]),3,2,3,iso=F,SigInv =cov(as.matrix(tt[,1:4])) )


###



rs=matrix(NA,4,9)

for(i in 1:4){
  g1=kmeans_inference_1f(as.matrix(t1),k = 3, 1,2,
                         feat = i, iso = F,
                         covMat=cov(t1), 
                         iter.max = 15)
  rs[i,1]=g1$pval
  
  g1=kmeans_inference_1f(as.matrix(t1),k = 3, 1,3,
                         feat = i, iso = F,
                         covMat=cov(t1), 
                         iter.max = 15)
  rs[i,2]=g1$pval
  
  g1=kmeans_inference_1f(as.matrix(t1),k = 3, 2,3,
                         feat = i, iso = F,
                         covMat=cov(t1), 
                         iter.max = 15)
  
  rs[i,3]=g1$pval
  
  g1=kmeans_inference_1f(as.matrix(tt[tt$ct==1,]),k = 2, 1,2,
                         feat = i, iso = F,
                         covMat=cov(tt[tt$ct==1,]), 
                         iter.max = 15)
  
  rs[i,4]=g1$pval
  
  g1=kmeans_inference_1f(as.matrix(tt[tt$ct==2,]),k = 2, 1,2,
                         feat = i, iso = F,
                         covMat=cov(tt[tt$ct==2,]), 
                         iter.max = 15)
  
  rs[i,5]=g1$pval
  
  g1=kmeans_inference_1f(as.matrix(tt[tt$ct==3,]),k = 2, 1,2,
                         feat = i, iso = F,
                         covMat=cov(tt[tt$ct==3,]), 
                         iter.max = 15)
  
  rs[i,6]=g1$pval
  
  g1=kmeans_inference_1f(as.matrix(tt),k = 3, 1,2,
                         feat = i, iso = F,
                         covMat=cov(tt), 
                         iter.max = 15)
  rs[i,7]=g1$pval
  
  g1=kmeans_inference_1f(as.matrix(tt),k = 3, 1,3,
                         feat = i, iso = F,
                         covMat=cov(tt), 
                         iter.max = 15)
  rs[i,8]=g1$pval
  g1=kmeans_inference_1f(as.matrix(tt),k = 3, 2,3,
                         feat = i, iso = F,
                         covMat=cov(tt), 
                         iter.max = 15)
  rs[i,9]=g1$pval
  
}
