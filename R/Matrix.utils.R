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
#' test<-data.frame(orderNum=as.factor(sample(1e6, 1e7, TRUE)),
#'    sku=as.factor(sample(1e3, 1e7, TRUE)),
#'    customer=as.factor(sample(1e4,1e7,TRUE)),
#'    state = sample(letters, 1e7, TRUE), amount=runif(1e7))
#' # For simple aggregations resulting in small tables, dcast.data.table (and even reshape2) will be faster
#' system.time(a<-dcast.data.table(as.data.table(test),sku~state,sum,value.var = 'amount')) # .5 seconds
#' system.time(b<-reshape2::dcast(test,sku~state,sum,value.var = 'amount')) # 2.61 seconds
#' system.time(c<-dcast.Matrix(test,sku~state,value.var = 'amount')) # 28 seconds
#' # However, this situation changes as the result set becomes larger
#' system.time(a<-dcast.data.table(as.data.table(test),customer~sku,sum,value.var = 'amount')) #4 .4 seconds
#' system.time(b<-reshape2::dcast(test,customer~sku,sum,value.var = 'amount')) # 34.7 seconds
#' system.time(c<-dcast.Matrix(test,customer~sku,value.var = 'amount')) # 27 seconds
#' # More complicated:
#' system.time(a<-dcast.data.table(as.data.table(test),customer~sku+state,sum,value.var = 'amount')) # 18.1 seconds, object size = 2084 Mb
#' system.time(b<-reshape2::dcast(test,customer~sku+state,sum,value.var = 'amount')) # Does not return
#' system.time(c<-dcast.Matrix(test,customer~sku:state,value.var = 'amount')) # 30.69 seconds, object size = 115.4 Mb
#' 
#' system.time(a<-dcast.data.table(as.data.table(test),orderNum~sku,sum,value.var = 'amount')) # Does not return
#' system.time(c<-dcast.Matrix(test,orderNum~sku,value.var = 'amount')) # 36.33 seconds, object size = 175Mb
#' 
#' system.time(c<-dcast.Matrix(test,orderNum~customer,value.var = 'amount')) # 31.28 seconds, object size = 176Mb
#' 
#' 
dcast.Matrix<-function(data,formula,value.var=NULL)
{
  alltms<-terms(formula,data=data)
  response<-all.vars(formula)[attr(alltms,"response")]
  tm<-attr(alltms,"term.labels")
  interactionsIndex<-grep(':',tm)
  interactions<-tm[interactionsIndex]
  simple<-setdiff(tm,interactions)
  i2<-strsplit(interactions,':')
  newterms<-unlist(lapply(i2,function (x) paste("paste(",paste(x,collapse=','),",","sep='_'",")")))
  newterms<-c(simple,newterms)
  newformula<-as.formula(paste('~0+',paste(newterms,collapse='+')))
  result<-sparse.model.matrix(newformula,data)
  values<-1
  if(!is.null(value.var))
    values<-data[,value.var]
  result<-result*values
  if(isTRUE(response>0))
  {
    data[,response]<-as.factor(data[,response])
    result<-aggregate.Matrix(result,data[,response,drop=FALSE])
  }
  return(result)
}

#' @examples
#' ## Not run:
#' system.time(d<-aggregate.Matrix(test[,'amount',drop=FALSE],test$orderNum))
aggregate.Matrix<-function(data,groupings,fun=colSums)
{
  if(!is(data,'Matrix'))
    data<-Matrix(as.matrix(data),sparse=TRUE)
  data<-Matrix(data)
  mapping<-sparse.model.matrix(~0+.,data.frame(groupings))
  result<-t(mapping) %*% data
}