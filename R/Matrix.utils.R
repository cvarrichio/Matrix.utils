#' Row-based functions for R objects
#'
#' Rowr allows the manipulation of R objects as if they were organized rows in a
#' way that is familiar to people used to working with databases.  It allows
#' more consistent and predictable output to common functions, and generalizes a
#' number of utility functions to to be failsafe with any number of objects.
#' @name Matrix.utils
#' @docType package
NULL


#' Robust alternative to \code{\link{Vectorize}} function that accepts any function with two 
#' or more arguments.  Returns a function that will work an arbitrary number of vectors, lists or 
#' data frames, though output may be unpredicatable in unusual applications.  The 
#' results are also intended to be more intuitive than Vectorize.
#' 
#' @param fun a two or more argument function
#' @param type like \code{MARGIN} in \code{\link{apply}}, except that \code{c(1,2)} is
#'   represented as a \code{3} instead.  By default, will \code{Reduce} single dimensional
#'   data handle everything else row-wise.
#' @export
#' @examples
#' ## Not run:
#' orders<-data.frame(orderNum=as.factor(sample(1e6, 1e7, TRUE)),
#'    sku=as.factor(sample(1e3, 1e7, TRUE)),
#'    customer=as.factor(sample(1e4,1e7,TRUE)),
#'    state = sample(letters, 1e7, TRUE), amount=runif(1e7))
#' # For simple aggregations resulting in small tables, dcast.data.table (and even reshape2) will be faster
#' system.time(a<-dcast.data.table(as.data.table(orders),sku~state,sum,value.var = 'amount')) # .5 seconds
#' system.time(b<-reshape2::dcast(orders,sku~state,sum,value.var = 'amount')) # 2.61 seconds
#' system.time(c<-dcast.Matrix(orders,sku~state,value.var = 'amount')) # 28 seconds
#' # However, this situation changes as the result set becomes larger
#' system.time(a<-dcast.data.table(as.data.table(orders),customer~sku,sum,value.var = 'amount')) #4 .4 seconds
#' system.time(b<-reshape2::dcast(orders,customer~sku,sum,value.var = 'amount')) # 34.7 seconds
#' system.time(c<-dcast.Matrix(orders,customer~sku,value.var = 'amount')) # 27 seconds
#' # More complicated:
#' system.time(a<-dcast.data.table(as.data.table(orders),customer~sku+state,sum,value.var = 'amount')) # 18.1 seconds, object size = 2084 Mb
#' system.time(b<-reshape2::dcast(orders,customer~sku+state,sum,value.var = 'amount')) # Does not return
#' system.time(c<-dcast.Matrix(orders,customer~sku:state,value.var = 'amount')) # 30.69 seconds, object size = 115.4 Mb
#' 
#' system.time(a<-dcast.data.table(as.data.table(orders),orderNum~sku,sum,value.var = 'amount')) # Does not return
#' system.time(c<-dcast.Matrix(orders,orderNum~sku,value.var = 'amount')) # 36.33 seconds, object size = 175Mb
#' 
#' system.time(c<-dcast.Matrix(orders,orderNum~customer,value.var = 'amount')) # 31.28 seconds, object size = 176Mb
#' 
#' 
dcast.Matrix<-function(data,formula,value.var=NULL,fun='sum')
{
  require(dplyr)
  #browser()
  alltms<-terms(formula,data=data)
  response<-rownames(attr(alltms,'factors'))[attr(alltms,'response')]
  tm<-attr(alltms,"term.labels")
  interactionsIndex<-grep(':',tm)
  interactions<-tm[interactionsIndex]
  simple<-setdiff(tm,interactions)
  i2<-strsplit(interactions,':')
  newterms<-unlist(lapply(i2,function (x) paste("paste(",paste(x,collapse=','),",","sep='_'",")")))
  newterms<-c(simple,newterms)
  newformula<-as.formula(paste(response,'~0+',paste(newterms,collapse='+')))
  result<-sparse.model.matrix(newformula,data,drop.unused.levels = TRUE)
  colnames(result)<-lapply(colnames(result),function (x) gsub('paste(',replacement='',x=x,fixed = TRUE) %>% gsub(pattern=', ',replacement='_',x=.,fixed=TRUE) %>% gsub(pattern='_sep = \"_\")',replacement='',x=.,fixed=TRUE))
  values<-1
  if(!is.null(value.var))
    values<-data[,value.var]
  result<-result*values
  if(isTRUE(response>0))
  {
    responses=all.vars(terms(as.formula(paste(response,'~0'))))
    result<-aggregate.Matrix(result,form=paste('~0 +',response),data[,responses,drop=FALSE],fun=fun)
  }
  return(result)
}

#' @examples
#' ## Not run:
#' system.time(d<-aggregate.Matrix(orders[,'amount',drop=FALSE],orders$orderNum))
aggregate.Matrix<-function(data,groupings=NULL,form=NULL,fun='sum')
{
  browser()
  if(!is(data,'Matrix'))
    data<-Matrix(as.matrix(data),sparse=TRUE)
  if(fun=='count')
    data<-data!=0
  groupings<-data.frame(groupings)
  groupings<-data.frame(lapply(groupings,as.factor))
  if(is.null(form))
    form<-as.formula('~0+.')
  form<-as.formula(form)
  mapping<-dcast.Matrix(groupings,form)
  result<-t(mapping) %*% data
  if(fun=='mean')
    result@x<-result@x/(aggregate.Matrix(data,groupings,fun='count'))@x
  return(result)
}

#' Joins two matrices
#' @examples
#' orders<-Matrix(as.matrix(data.frame(orderNum=1:1000,customer=sample(100,1000,TRUE))),sparse=TRUE)
#' cancelledOrders<-Matrix(as.matrix(data.frame(orderNum=sample(1000,100),cancelled='Y)),sparse=TRUE)
#' skus<-Matrix(as.matrix(data.frame(orderNum=sample(1000,10000,TRUE),sku=sample(1000,10000,TRUE),amount=runif(10000))),sparse=TRUE)
join.Matrix<-function(x,y,rownames.x=rownames(x),rownames.y=rownames(y))
{
  browser()
#   rownames.x<-as.character(rownames.x)
#   rownames.y<-as.character(rownames.y)
#   rownames(x)<-rownames.x
#   rownames(y)<-rownames.y
  #inter<-intersect(rownames.x,rownames.y)
  #y<-y[inter,,drop=FALSE]
  indices1<-match(rownames.x,rownames.y)
  indices1[is.na(indices1)]<-nrow(y)+1
  y2<-rbind2(y,NA)[indices1,]
  colnames(y2)[colnames(y2) %in% colnames(x)]<-paste('y',colnames(y2)[colnames(y2) %in% colnames(x)])
  result<-cbind2(x,y2)
  #indices2<-match(rownames.y,rownames.x)
  #result<-cbind2(x,y[rownames.x,])
  return(result)
}
