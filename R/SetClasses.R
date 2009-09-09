setClass("flowObj",contains="flowClust",representation=list(DATA="environment"),prototype=list(DATA=new.env(hash=TRUE,parent=emptyenv())),package="flowMerge");

setClass("flowMerge",contains="flowObj",representation=list(merged="numeric",entropy="numeric"),prototype=list(merged=0,entropy=Inf),package="flowMerge")

