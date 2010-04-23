map<-function(z,...){
apply(z,1,function(x){w<-(which(x==max(x)));if(length(w)!=0){w}else{NA}})
}
setMethod("merge",signature=signature(x="flowObj",y="missing"),function(x,...){
  k <- x@K;
  resultObject <- list()
  resultObject[[k]] <- as(x,"flowMerge");
  resultObject[[k]]@entropy <- -2*sum(resultObject[[k]]@z*log(resultObject[[k]]@z,base=2), na.rm=T)
  d <- dim(resultObject[[k]]@mu)[1];
  s <- apply(resultObject[[k]]@sigma,1,function(x){xx<-eigen(x);list(solve(xx$vectors%*%sqrt(diag(xx$values))%*%t(xx$vectors)));});
  #resultObject[[k]]@ssd<- 2*sum(sapply(1:d,function(i)sapply(1:d,function(j)sqrt(sum((s[[i]][[1]]%*%resultObject[[k]]@mu[i,]-s[[j]][[1]]%*%resultObject[[k]]@mu[j,])^2)))))

  if(k > 2){
    for (kk in (k-1):2) resultObject[[kk]] <- mergeClusters(resultObject[[kk+1]], getData(resultObject[[k]]))
    resultObject[[1]] <- mergeClusters(resultObject[[2]], getData(resultObject[[k]]))
  }else if(k==2){
    resultObject[[1]] <- mergeClusters(resultObject[[2]], getData(resultObject[[k]])) 
  }
    resultObject <- lapply(resultObject, updateU);
    resultObject <- lapply(resultObject, flagOutliers);
    resultObject<-lapply(resultObject,function(x){ruleOutliers(x)<-list(level=0.9);x});
  resultObject;
})

setMethod("split",signature=signature(f="flowMerge",x="missing"),function(f, drop = FALSE, population = NULL, split = NULL, rm.outliers = TRUE, ...){
    split(f = as(f, "flowClust"), x = getData(f), population = population, split = split, rm.outliers = rm.outliers, ...);
});

setGeneric("updateU", function(object){standardGeneric("updateU")});
setMethod("updateU", signature = signature(object = "flowMerge"), function(object){
  p<-object@K;
  q <- object@varNames;
  #q <- dim(o@mu)[2];
  if(length(object@lambda)==1){#test whether the transformation parameter is 1  (length 0)
    dprime <- flowClust::box(exprs(object@DATA[["1"]])[,q],object@lambda);
  }else{
    dprime<-exprs(object@DATA[["1"]])[,q]
  }   
  u<-sapply(1:object@K,function(i){s<-solve(object@sigma[i,,]);apply(dprime,1,function(x)(x-object@mu[i,])%*%s%*%(x-object@mu[i,]))})
  if(!is.infinite(object@nu)){
      object@u <- (length(q)+object@nu)/(u+object@nu)
  }else{
    object@u<-u;
  }
  object@label <- map(object@z);
  object;
})

setMethod("getData",signature=signature(obj="flowObj"),function(obj){obj@DATA[["1"]]});
setMethod("getData",signature=signature(obj="flowMerge"),function(obj){getData(as(obj,"flowObj"))});
setMethod("plot",signature=signature(x="flowObj",y="missing"),
          function(x,new.window=TRUE,...){
            #TODO include some code so that the dimension names are chosen correctly.
            comb<-as.matrix(combn(length(x@varNames),2));
            if(new.window){
            par(mfrow=c(ceiling(sqrt(dim(comb)[2])),ceiling(sqrt(dim(comb)[2]))))
            }
            if(length(x@varNames)>2){
                apply(comb,2,function(i){
                    selectMethod("plot",signature=c(x="flowClust",y="missing"))(as(x,"flowClust"),data=getData(x),subset=i,xlab=parameters(getData(x))$desc[i[1]+2],ylab=parameters(getData(x))$desc[i[2]+2],...);
             })
            }else{
                selectMethod("plot",signature=c(x="flowClust",y="missing"))(as(x,"flowClust"),data=getData(x),xlab=parameters(getData(x))$desc[comb[1]],ylab=parameters(getData(x))$desc[comb[2]],...);
            }
          })

setMethod("plot",signature=signature(x="flowMerge",y="missing"),
          function(x,new.window=FALSE,...){
            #TODO include some code so that the dimension names are chosen correctly.
            comb<-as.matrix(combn(length(x@varNames),2));
            if(new.window){
            par(mfrow=c(ceiling(sqrt(dim(comb)[2])),ceiling(sqrt(dim(comb)[2]))))
            }
            if(length(x@varNames)>2){
                apply(comb,2,function(i){
                    selectMethod("plot",signature=c(x="flowClust",y="missing"))(as(x,"flowClust"),data=getData(x),subset=i,xlab=parameters(getData(x))$desc[i[1]],ylab=parameters(getData(x))$desc[i[2]],...);
             })
            }else{
                selectMethod("plot",signature=c(x="flowClust",y="missing"))(as(x,"flowClust"),data=getData(x),xlab=parameters(getData(x))$desc[comb[1]],ylab=parameters(getData(x))$desc[comb[2]],...);
            }
          })
          
setGeneric("fitPiecewiseLinreg",function(x,plot=FALSE,...){standardGeneric("fitPiecewiseLinreg")})

setMethod("fitPiecewiseLinreg",signature=signature(x="list"),{
    function(x,plot=FALSE,...){
        if(!all(unlist(lapply(x,function(y)is(y,"flowMerge"))))){
            stop("x is not a valid list of flowMerge objects");        
        }
        entropy<-flowMerge:::ENT(x);
        N<-(unlist(lapply(x,function(y)y@merged)))
        N<-cumsum(c(0,N[-length(N)]))
        l<-length(entropy);
    
        #Two cases.. if there's more than two clusters, or if there's just two clusters. 
        #If there are more than two clusters.. we can fit a changepoint model.
        if(l>3){
            #Try all positions for the changepoint, from 2..(l-1), and compute the sum of squared residuals
            r<-sapply(2:(l-1),function(b) try(sum(lm(entropy~N,subset=c(1:b))$residuals^2)+sum(lm(entropy~N,subset=c(b:l))$residuals^2)))
            r<-as.numeric(r);
            r2<-sum(lm(entropy~N)$residuals^2);
            bic<-c(l*log(r2/l)+2*log(l),l*log(r/l)+5*log(l));        
            m<-which.min(bic);        
            if(plot){
                c1<-coefficients(lm(entropy~N,subset=c(1:m)));c2<-coefficients(lm(entropy~N,subset=c(m:l)))
            color<-rep(1,l);
            color[m]<-2;
            plot(N,entropy,col=color,pch=20,main="Entropy of Clustering",xlab="Cumulative Number of Merged Observations",ylab="Entropy");
            lines(N[1:m],N[1:m]*c1[2]+c1[1],col="red");
            lines(N[m:l],N[m:l]*c2[2]+c2[1],col="red");         
            }
            return(m);
            if(m==1){
                warning("No changepoint found. Returning max BIC solution");
                return(l);
            }else{
                return(m);
            }        
        }else if(l==2){
            warning("Two clusters.. no changepoint possible. Returning max BIC solution");
            return(2);
        }else if(l==3)
        {
            r2<-sum(lm(entropy~N)$residuals^2);
            r<-sapply(2:(l-1),function(b) try(list(lm(entropy~N,subset=c(1:b)), lm(entropy~I(1:l),subset=c(b:l)))));
            a<-coefficients(r[[1]])[2];
            b<-coefficients(r[[2]])[2];
            angle<-atan(abs(a-b)/(1+a*b))*180/pi
            if(angle>1){
                warning("Possible changepoint detected by angle between line segments")
                return(2)                            
            }else{
                warning("No changepoint detected, returning max bic solution")
                return(3);            
            }
        }else if(l==1){
         return(1);   
        }
    }
})   

setMethod("summary",signature=signature("flowObj"),
function(object,...){
    summary(as(object,"flowClust"));
    print(summary(object@uncertainty));
    cat("Data summary:\n");
    summary(getData(object));    
})
setMethod("summary",signature=signature("flowMerge"),
function(object,...){
    cat("** Merging Information **\n");
    cat("Entropy of clustering: ",object@entropy,"\n");
    cat("Normalized Entropy of clustering: ",NENT(list(object)),"\n");
    summary(as(object,"flowObj"));
})

setGeneric("flagOutliers",function(object,...){standardGeneric("flagOutliers")});
setMethod("flagOutliers",signature=signature("flowMerge"),
function(object,...){
    if(object@ruleOutliers[1]==0){        
        object@flagOutliers<-sapply(1:nrow(object@u),function(i){
        object@u[i,object@label[i]]<object@ruleOutliers[2]||object@z[i,object@label[i]]<object@ruleOutliers[3]})
    }else if(object@ruleOutliers[1]==1){
        object@flagOutliers<-sapply(1:nrow(object@u),function(i){
        object@flagOutliers<-object@z[i,object@label[i]]<object@ruleOutliers[3]})
    }
    object;
})

setMethod("show",signature=signature("flowObj"),
function(object){
       cat("Object of class 'flowObj'\n");
       cat("This object has the following slots:\n");
       cat("DATA, expName, varNames, K, w, mu, sigma, lambda, nu, z, u, label, uncertainty, ruleOutliers, flagOutliers, rm.min, rm.max, logLike, BIC, ICL\n");
});

setMethod("show",signature=signature("flowMerge"),
function(object){
       cat("Object of class 'flowMerge'\n");
       cat("This object has the following slots:\n");
       cat("merged, entropy, DATA, expName, varNames, K, w, mu, sigma, lambda, nu, z, u, label, uncertainty, ruleOutliers, flagOutliers, rm.min, rm.max, logLike, BIC, ICL\n");
});
