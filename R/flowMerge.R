
checkForRemoteErrors <- function(val)
    {
        count <- 0
    firstmsg <- NULL
    for (v in val) {
        if (inherits(v, "try-error")) {
            count <- count + 1
            if (count == 1) 
                firstmsg <- v
        }
    }
    if (count == 1) 
        print(paste("one node produced an error: ", firstmsg))
    else if (count > 1) 
        print(paste(count, " nodes produced errors; first error: ", 
            firstmsg))
    val
}

initPFlowMerge<-function(cl){
    if(!is(cl,"cluster")){
        stop("cl must be a SNOW cluster");
    }
    r<-clusterEvalQ(cl,library(flowMerge));
    if(any(unlist(lapply(r,function(x)is(x,"try-error"))))){
        print(r);
     stop("Some cluster nodes failed to initialize correctly");   
    }
}

pFlowClust<-function(flowData,cl,K=1:15,B.init=100,tol.init=1e-2,tol=1e-5,B=1000,randomStart=50,nu=4,nu.est=1,trans=1,varNames=NA){
if(any(is.na(varNames))){
stop("varNames must be defined")
}
if(!is(cl,"cluster")){
stop("cl must be an object of type cluster from the snow package");
}
if(is(flowData,"list")){
if(all(unlist(lapply(flowData,function(x)is(x,"flowFrame"))))){
result<-lapply(flowData,function(d)clusterMap(cl,function(...)try(flowClust(...)),list(d),varNames=list(varNames),K=K,B=B,tol.init=tol.init,tol=tol,B.init=B.init,randomStart=randomStart,nu=nu,nu.est=nu.est,trans=trans));
}}
else if(is(flowData,"flowSet")){
result<-fsApply(flowData,function(d)clusterMap(cl,function(...)try(flowClust(...)),list(d),varNames=list(varNames),K=K,B=B,tol.init=tol.init,tol=tol,B.init=B.init,randomStart=randomStart,nu=nu,nu.est=nu.est,trans=trans));
}
else if(is(flowData,"flowFrame")){
result<-clusterMap(cl,function(...)try(flowClust(...)),list(flowData),varNames=list(varNames),K=K,B=B,tol.init=tol.init,tol=tol,B.init=B.init,randomStart=randomStart,nu=nu,nu.est=nu.est,trans=trans)
}
else{
stop("flowData should be a flowFrame, list of flowFrames, or flowSet");
}
result;
}

pFlowMerge<-function(flowData,cl,K=1:15,B.init=100,tol.init=1e-2,tol=1e-5,B=500,randomStart=10,nu=4,nu.est=0,trans=1,varNames=NA){
    if(!is.null(cl)){
        initPFlowMerge(cl);
        result<-pFlowClust(flowData,cl,varNames=varNames,K=K,B.init=B.init,tol.init=tol.init,tol=tol,B=B,randomStart=randomStart,nu=nu,nu.est=nu.est,trans=trans);    
    }
    if(is.null(cl)){
        if(is(flowData,"list")|is(flowData,"flowSet")){
            result<-lapply(as(flowData,"list"),function(x)mapply(function(...)try(flowClust(...)),list(x),varNames=list(varNames),K=K,B.init=B.init,tol.init=tol.init,tol=tol,B=B,randomStart=randomStart,nu=nu,nu.est=nu.est,trans=trans));
        }else if(is(flowData,"flowFrame")){
          result<-mapply(function(...)try(flowClust(...)),list(flowData),varNames=list(varNames),K=K,B.init=B.init,tol.init=tol.init,tol=tol,B=B,randomStart=randomStart,nu=nu,nu.est=nu.est,trans=trans)
        }else{
            stop("flowData must be a flowFrame, list of flowFrames or flowSet")
        }   
    }
    ##Extract the max BIC objects
    if(is(result,"list")){##This nested code tests if the output is a list of lists of flowClust results. Inelegant, maybe, but that's R for you.
        if(all(unlist(lapply(result,function(x)is(x,"list"))))){
            if(all(unlist(lapply(result,function(x)lapply(x,function(y)is(y,"flowClust")))))){
                tmp<-lapply(result,function(x)x[[which.max(BIC(x))]]);
                ##result should be a list of flowClust objects, one per flowFrame input
                o<-sapply(1:length(tmp),function(i)flowObj(tmp[[i]],flowData[[i]]));
                m<-lapply(o,function(x)merge(x));
                result<-lapply(m,function(x){i<-fitPiecewiseLinreg(x,plot=T);x[[i]]});
            }else{
                stop("The result must be the output from a parallel call to flowClust");   
            }
        }else{
            if(all(unlist(lapply(result,function(x)is(x,"flowClust")|is(x,"try-error"))))){
                o<-flowObj(result[[which.max(BIC(result))]],flowData);
                m<-merge(o);
                i<-fitPiecewiseLinreg(m,plot=T);
                result<-m[[i]];
            }else{
                stop("oops! The input to pFlowMerge from pFlowClust is incorrect.. neither a list of lists of flowClust objects, nor a list of flowClust objects");   
            }   
        }
    }    
    result;
}


#cl<-makeSOCKcluster(c(rep("finakg@got04",4),rep("recherche@gotsrv02",8),rep("recherche@gotsrv03",8)))
#clusterEvalQ(cl,library(flowClust))
#clusterEvalQ(cl,source("/Users/recherche/merge.R"));

flowObj<-function(flowC=NULL,flowF=NULL){
  if(!is(flowC,"flowClust")){
    stop("flowC must be a flowClust object");
  }
  if(!is(flowF,"flowFrame")){
    stop("flowF must be a flowFrame object");
  }
  o <- new("flowObj",flowC);
  o@DATA=new.env(hash=TRUE,parent=emptyenv());
  o@DATA[["1"]]<-flowF;
  o;
}

mergeClusters2 <- function(object, a, b){
#    ly <- nrow(object@z)
    #py <- ncol(data)
	py<-dim(object@mu)[1]
	dy<-dim(object@mu)[2]
    
    # update mu, sigma, w
    P <- object@w[a]+object@w[b];
    MU<-(object@w[a]*object@mu[a,]+object@w[b]*object@mu[b,])/P
    if(is.infinite(object@nu)){
        SIGMA<-(object@w[a]*(object@mu[a,]%*%t(object@mu[a,])+(object@sigma[a,,])+object@w[b]*(object@mu[b,]%*%t(object@mu[b,])+(object@sigma[b,,]))/(P) - MU%*%t(MU)))
    }else{
    SIGMA<-((object@w[a]*(object@mu[a,]%*%t(object@mu[a,])+((object@nu)/(object@nu-2))*object@sigma[a,,])+object@w[b]*(object@mu[b,]%*%t(object@mu[b,])+((object@nu)/(object@nu-2))*object@sigma[b,,]))/(P) - MU%*%t(MU))*((object@nu-2)/object@nu)
    }
    object@w[a]<-P; object@w<-object@w[-b];
    object@mu[a,]<-MU;
    dims<-dim(object@mu);
    object@z[,a] <- object@z[,a] + object@z[,b]
    object@sigma[a,,]<-SIGMA;
	sigma<-object@sigma[-b,,];
	if(!is.array(sigma)&dy==1){
		#One dimensional case.. need to make this a proper array
		object@sigma<-array(data=sigma,c(py-1,dy,dy))
	}
    #if sigma is a single covariance matrix, then set it up into an appropriate array
    if(length(dim(object@sigma))==2){
        d<-dim(object@sigma);
        object@sigma<-array(object@sigma,dim=c(1,d[1],d[2]));
    }
    if(dims[1]==2){
      object@mu <- t(matrix(object@mu[-b,]));
	  object@z <- matrix(object@z[,-b]);
    }else if(dims[2]==1){
		object@mu<-matrix(object@mu[-b,])
		object@z <- object@z[,-b];
	}else {
      object@mu<-object@mu[-b,];
      object@z <- object@z[,-b];
    }
  
    


    # update K, ICL
    object@K <- object@K-1
    
    object@entropy <- -2*sum(object@z*log(object@z,base=2), na.rm=T)  
    object@ICL <- object@BIC + 2*sum(object@z*log(object@z),na.rm=T)
    return(object)
}

.computeDeltaE<-function(object,a,b){
    list(a=a,b=b,e=-2*sum(object@z[,c(a,b)]*log(object@z[,c(a,b)]),na.rm=TRUE));
}

.compMD<-function(object,a,b){
    list(a=a,b=b,e=mahalanobis(object@mu[a,],object@mu[b,],object@sigma[b,,])+mahalanobis(object@mu[b,],object@mu[a,],object@sigma[a,,]))
}

.combFunc1<-function(x,y){
    list(x,y)[[which.max(c(x$e,y$e))]]
}
.combFunc2<-function(x,y){
    list(x,y)[[which.min(c(x$e,y$e))]]
}

mergeClusters <- function(object,metric="entropy") {
    K <- object@K
    minEnt <- Inf
    maxPLL<- -Inf
	a<-0;
	b<-0;
    if(!is.na(match("doMC",installed.packages()))){
		require(doMC)
        doMC::registerDoMC();
    }else if(!is.na(match("doSNOW",installed.packages()))){
		require(doSNOW)
        doSNOW::registerDoSNOW()
    }
    if(any(grepl("doSNOW",loadedNamespaces()))|any(grepl("doMC",loadedNamespaces()))){
    	if(metric=="entropy"){
    		resObject<-foreach (a = 1:(K-1),.combine=.combFunc1, .options.multicore=list(preschedule=TRUE))%dopar%{ 
     			foreach (b = (a+1):K,.combine=.combFunc1,.options.multicore = list(preschedule=TRUE)) %dopar%{
      				.computeDeltaE(object,a,b);
				}
        
    		}
		}
		else if(metric=="mahalanobis"){
        	resObject<-foreach (a = 1:(K-1),.combine=.combFunc2, .options.multicore=list(preschedule=TRUE))%dopar%{ 
     			foreach (b = (a+1):K,.combine=.combFunc2,.options.multicore = list(preschedule=TRUE)) %dopar%{
         			.compMD(object,a,b);
        		}
      		}
    	}    
		resObject
	}
	else{
    	if(metric=="entropy"){
    		resObject<-foreach (a = 1:(K-1),.combine=.combFunc1, .options.multicore=list(preschedule=TRUE))%do%{ 
     			foreach (b = (a+1):K,.combine=.combFunc1,.options.multicore = list(preschedule=TRUE)) %do%{
        			.computeDeltaE(object,a,b);
      			}
    		}
		}
		else if(metric=="mahalanobis"){
        	resObject<-foreach (a = 1:(K-1),.combine=.combFunc2, .options.multicore=list(preschedule=TRUE))%do%{ 
     			foreach (b = (a+1):K,.combine=.combFunc2,.options.multicore = list(preschedule=TRUE)) %do%{
         			.compMD(object,a,b);
        		}
      		}
    	}    
	}
}

BIC <- function(x)as.numeric(unlist(lapply(x,function(q)try(q@BIC,silent=TRUE))))
ICL <- function(x)as.numeric(unlist(lapply(x,function(q)try(q@ICL,silent=TRUE))))
ENT <-function(x)as.numeric(unlist(lapply(x,function(q)try(q@entropy,silent=TRUE))))
NENT<-function(x)as.numeric(unlist(lapply(x,function(q)try(q@entropy/(q@K*dim(q@z)[1]),silent=TRUE))))




