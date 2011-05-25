setClass("flowObj",contains="flowClust",representation=list(DATA="environment"),prototype=list(DATA=new.env(hash=TRUE,parent=emptyenv())),package="flowMerge");

setClass("flowMerge",contains="flowObj",representation=list(merged="numeric",entropy="numeric",mtree="graphNEL"),prototype=list(merged=0,entropy=Inf,mtree=new("graphNEL",edgemode="directed")),package="flowMerge")

