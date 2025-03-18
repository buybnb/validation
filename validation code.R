### Function call validation

##validation(discovery,validation)
###d is the discovery data which only including features you plan to do the clustering
###v is the validation data which has exactly same features
### k is number of clusters


validation=function(d,v,k){
  
  zx=combn(k,2)
  
  rsa = matrix(NA,2*ncol(zx)+2*k,3)
  rsa[,3]="-"
  res_Kmeans = KmeansInference::kmeans_estimation(as.matrix(d), k, seed = NULL)
  cl1 = res_Kmeans$final_cluster
  
  ##Di
  
  for(i in 1:k){
  
    rsa[i,1] = paste("indivisible_D",i,sep="")
    
    rsa[i,2]==kmeans_inference(d[cl1==i,],2,1,2,iso=F,SigInv =solve(cov(d[cl1==i,])),seed=NULL)$pval

  }

  for(i in 1:ncol(zx)){
  
    rsa[k+i,1]=paste("nonmergeable_D",zx[1,i]," vs D",zx[2,i],sep="")

    rsa[k+i,2]  =kmeans_inference(d,k,zx[1,i],zx[2,i],iso=F,SigInv =solve(cov(d)),seed=NULL)$pval

  }

  zz=clus(res_Kmeans,v)
  cl2=zz$val
  
  xx=rbind(d,v)
  ct = c(cl1,cl2)
  
  for(i in 1:k){
  
    rsa[k+ncol(zx)+i,1] = paste("indivisible_D",i," U V",i,sep="")
  
    rsa[k+ncol(zx)+i,2]=kmeans_inference(xx[ct==i,],2,1,2,iso=F,SigInv =solve(cov(xx[ct==i,])),seed=NULL)$pval

    rsa[k+ncol(zx)+i,3]=ifelse(rsa[k+ncol(zx)+i,2]<0.05,"No","Yes")
  
  }
  
  for(i in 1:ncol(zx)){
    
    rsa[2*k+ncol(zx)+i,1]=paste("nonmergeable_D",zx[1,i]," U V",zx[1,i]," vs D",zx[2,i]," U V",zx[2,i],sep="")
    
    rsa[2*k+ncol(zx)+i,2]=kmeans_inference(xx[ct%in%zx[,i],],2,1,2,iso=F,SigInv =solve(cov(xx[ct%in%zx[,i],])),seed=NULL)$pval
    
  }
  
  colnames(rsa)=c("Criteria","p-value","if validate")
  rsa
}


#validation(d,v,3)
