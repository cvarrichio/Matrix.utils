#' Data.frame-like operations on sparse and dense \code{Matrix} objects,
#' 
#' Implements cast, aggregate, and join for \code{\link{Matrix}} and matrix-like
#' objects.
#' 
#' @name Matrix.utils
#' @docType package
#' @import Matrix
#' @importFrom stats aggregate
NULL

#' Casts or pivots a long \code{data frame} into a wide sparse matrix.
#' 
#' Similar in function to \code{\link[reshape2]{dcast}}, but produces a sparse 
#' \code{\link{Matrix}} as an output. Sparse matrices are beneficial for this 
#' application because such outputs are often very wide and sparse. Conceptually
#' similar to a \code{pivot} operation.
#' 
#' Casting formulas are slightly different than those in \code{dcast} and follow
#' the conventions of \code{\link{model.matrix}}. See \code{\link{formula}} for 
#' details.  Briefly, the left hand side of the \code{~} will be used as the 
#' grouping criteria.  This can either be a single variable, or a group of 
#' variables linked using \code{:}.  The right hand side specifies what the 
#' columns will be. Unlike \code{dcast}, using the \code{+} operator will append
#' the values for each variable as additional columns.  This is useful for 
#' things such as one-hot encoding.  Using \code{:} will combine the columns as 
#' interactions.
#' 
#' @param data a data frame
#' @param formula casting \code{\link[stats]{formula}}, see details for
#'   specifics.
#' @param fun.aggregate name of aggregation function.  Defaults to 'sum'
#' @param value.var name of column that stores values to be aggregated
#' @seealso \code{\link[reshape]{cast}}
#' @seealso \code{\link[reshape2]{dcast}}
#' @seealso \code{\link[data.table]{dcast.data.table}}
#' @export
#' @examples
#' 
#' #Classic air quality example
#' melt<-function(data,idColumns)
#' {
#'   cols<-setdiff(colnames(data),idColumns)
#'   results<-lapply(cols,function (x) cbind(data[,idColumns],variable=x,value=as.numeric(data[,x])))
#'   results<-Reduce(rbind,results)
#' }
#' names(airquality) <- tolower(names(airquality))
#' aqm <- melt(airquality, idColumns=c("month", "day"))
#' dMcast(aqm, day:month ~variable,fun.aggregate = 'mean',value.var='value')
#' dMcast(aqm, month ~ variable, fun.aggregate = 'mean',value.var='value') 
#' 
#' #One hot encoding
#' #Preserving numerics
#' dMcast(warpbreaks,~.)
#' #Pivoting numerics as well
#' 
#' \dontrun{
#' orders<-data.frame(orderNum=as.factor(sample(1e6, 1e7, TRUE)), 
#'    sku=as.factor(cast.Matrsample(1e3, 1e7, TRUE)), 
#'    customer=as.factor(sample(1e4,1e7,TRUE)), 
#'    state = sample(letters, 1e7, TRUE),
#'    amount=runif(1e7)) 
#' # For simple aggregations resulting in small tables, dcast.data.table (and
#' # even reshape2) will be faster
#' system.time(a<-dcast.data.table(as.data.table(orders),sku~state,sum,
#'    value.var = 'amount')) # .5 seconds 
#' system.time(b<-reshape2::dcast(orders,sku~state,sum,
#'    value.var = 'amount')) # 2.61 seconds 
#' system.time(c<-dMcast(orders,sku~state,
#'    value.var = 'amount')) # 28 seconds 
#'    
#' # However, this situation changes as the result set becomes larger 
#' system.time(a<-dcast.data.table(as.data.table(orders),customer~sku,sum,
#'    value.var = 'amount')) # 4.4 seconds 
#' system.time(b<-reshape2::dcast(orders,customer~sku,sum,
#'    value.var = 'amount')) # 34.7 seconds 
#'  system.time(c<-dMcast(orders,customer~sku,
#'    value.var = 'amount')) # 27 seconds 
#'    
#' # More complicated: 
#' system.time(a<-dcast.data.table(as.data.table(orders),customer~sku+state,sum,
#'    value.var = 'amount')) # 18.1 seconds, object size = 2084 Mb 
#' system.time(b<-reshape2::dcast(orders,customer~sku+state,sum,
#'    value.var = 'amount')) # Does not return 
#' system.time(c<-dMcast(orders,customer~sku:state,
#'    value.var = 'amount')) # 30.69 seconds, object size = 115.4 Mb
#' 
#' system.time(a<-dcast.data.table(as.data.table(orders),orderNum~sku,sum,
#'    value.var = 'amount')) # Does not return 
#' system.time(c<-dMcast(orders,orderNum~sku,
#'    value.var = 'amount')) # 36.33 seconds, object size = 175Mb
#' 
dMcast<-function(data,formula,fun.aggregate='sum',value.var=NULL)
{
  browser()
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
  brokenNames<-grep('paste(',colnames(result),fixed = TRUE)
  colnames(result)[brokenNames]<-lapply(colnames(result)[brokenNames],function (x) {
    x<-gsub('paste(',replacement='',x=x,fixed = TRUE) 
    x<-gsub(pattern=', ',replacement='_',x=.,fixed=TRUE) 
    x<-gsub(pattern='_sep = \"_\")',replacement='',x=.,fixed=TRUE)
    return(x)
  })
  values<-1
  if(!is.null(value.var))
    values<-data[,value.var]
  result<-result*values
  if(isTRUE(response>0))
  {
    responses=all.vars(terms(as.formula(paste(response,'~0'))))
    result<-aggregate.Matrix(result,form=paste('~0 +',response),data[,responses,drop=FALSE],fun=fun.aggregate)
  }
  return(result)
}

#' Compute summary statistics of a Matrix.
#' 
#' Similar to \code{\link[stats]{aggregate}}.  Splits the matrix into groups as
#' specified by groupings, which can be one or more variables. Aggregation
#' function will be applied to all columns in data, or as specified in formula. 
#' Warning: groupings will be made dense if it is sparse, though data will not.
#' 
#' @param x a \code{\link{Matrix}} or matrix-like object
#' @param groupings an object coercible to a group of factors defining the groups
#' @param form \code{\link[stats]{formula}}
#' @param fun name of aggregation function to be applied to all columns in data
#' @seealso \code{\link[dplyr]{summarise}}
#' @seealso \code{\link[plyr]{summarise}}
#' @seealso \code{\link[stats]{aggregate}}
#' @export
#' @export aggregate.Matrix
#' @examples
#' skus<-Matrix(as.matrix(data.frame(
#'    orderNum=sample(1000,10000,TRUE),
#'    sku=sample(1000,10000,TRUE),
#'    amount=runif(10000))),sparse=TRUE)
#' a<-aggregate.Matrix(skus[,'amount'],skus[,'sku',drop=FALSE])
#' 
#' m<-rsparsematrix(1000000,100,.001)
#' labels<-as.factor(sample(1e4,1e6,TRUE))
#' b<-aggregate.Matrix(m,labels)
#' 
#' \dontrun{
#' orders<-data.frame(orderNum=as.factor(sample(1e6, 1e7, TRUE)),
#'    sku=as.factor(sample(1e3, 1e7, TRUE)),
#'    customer=as.factor(sample(1e4,1e7,TRUE)),
#'    state = sample(letters, 1e7, TRUE), amount=runif(1e7))
#' system.time(d<-aggregate.Matrix(orders[,'amount',drop=FALSE],orders$orderNum))
#' }
aggregate.Matrix<-function(x,groupings=NULL,form=NULL,fun='sum')
{
#   if(!is(x,'Matrix'))
#     x<-Matrix(as.matrix(x),sparse=TRUE)
  if(fun=='count')
    x<-x!=0
  if(is(groupings,'Matrix'))
    groupings<-as.matrix(groupings)
  groupings<-data.frame(groupings)
  groupings<-data.frame(lapply(groupings,as.factor))
  if(is.null(form))
    form<-as.formula('~0+.')
  form<-as.formula(form)
  mapping<-dMcast(groupings,form)
  result<-t(mapping) %*% x
  if(fun=='mean')
    result@x<-result@x/(aggregate.Matrix(x,groupings,fun='count'))@x
  return(result)
}

#'Merges two Matrices or matrix-like objects 
#'
#'Similar to \code{\link{merge}} or \code{\link[plyr]{join}. Rudimentary at this point; can do left and 
#'inner joins that only find the first matching row; equivalent to 
#'plyr::join(match='first').
#'
#'@param x \code{\link{Matrix}} or matrix-like object
#'@param y \code{Matrix} or matrix-like object
#'@param by.x vector indicating the names to match from \code{Matrix} x
#'@param by.y vector indicating the names to match from \code{Matrix} y
#'@param type type of join: currently on left and inner are supported
#'@export
#' @examples
#' orders<-Matrix(as.matrix(data.frame(orderNum=1:1000,
#'    customer=sample(100,1000,TRUE))),sparse=TRUE)
#' cancelledOrders<-Matrix(as.matrix(data.frame(orderNum=sample(1000,100),
#'    cancelled=1)),sparse=TRUE)
#' skus<-Matrix(as.matrix(data.frame(orderNum=sample(1000,10000,TRUE),
#'    sku=sample(1000,10000,TRUE),
#'    amount=runif(10000))),sparse=TRUE)
#' a<-join.Matrix(orders,cancelledOrders,orders[,'orderNum'],cancelledOrders[,'orderNum'])
#' b<-join.Matrix(orders,cancelledOrders,orders[,'orderNum'],cancelledOrders[,'orderNum'],type='inner')
#' c<-join.Matrix(orders,skus,orders[,'orderNum'],skus[,'orderNum'])
#' 
#' \dontrun{
#' ## Not run:
#' orders<-data.frame(orderNum=sample(1e6, 1e7, TRUE),
#'    sku=sample(1e3, 1e7, TRUE),
#'    customer=sample(1e4,1e7,TRUE))
#' cancelledOrders<-data.frame(data.frame(orderNum=sample(1e6,1e5),cancelled=1))
#' system.time(b<-join.Matrix(orders,cancelledOrders,orders[,'orderNum'],cancelledOrders[,'orderNum'],type='inner'))
#' #The following is the equivalent call in plyr, but returns an error due to a bug in plyr
#' system.time(c<-plyr::join(orders,cancelledOrders,type='inner',match='first')) 
#'}
join.Matrix<-function(x,y,by.x=rownames(x),by.y=rownames(y),type='left')
{
  indices1<-match(by.x,by.y)
  indices1[is.na(indices1)]<-nrow(y)+1
  y2<-rbind2(y,NA)[indices1,]
  colnames(y2)[colnames(y2) %in% colnames(x)]<-paste('y',colnames(y2)[colnames(y2) %in% colnames(x)])
  result<-cbind2(x,y2)
  if(type=='inner')
    result<-result[!(indices1==(nrow(y)+1)),]
  return(result)
}

melt<-function(data,idColumns)
{
  cols<-setdiff(colnames(data),idColumns)
  results<-lapply(cols,function (x) cbind(data[,idColumns],variable=x,value=as.numeric(data[,x])))
  results<-Reduce(rbind,results)
}