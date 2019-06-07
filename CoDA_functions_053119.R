#The compositional data analysis functions required to reproduce the analyses in Brown et al. 2019 mBio

#error-checking functions
format.check <- function(x){
  if (any(x==0)){
    stop('Zeroes detected, please replace')
  }
  if (is.vector(x)){
    m <- matrix(x, nrow = 1)
    colnames(m) <- names(x)
    m
  } else {
    x
  }
}

#geometric mean of rows, weighted or not
geomeans.r <- function(df, p=NULL, na.rm=FALSE){
  if (is.null(p)){
    exp(rowMeans(log(df), na.rm=na.rm))
  } else {
    exp((rowSums(log(df) %*% diag(p)))/sum(p))
  }
}

#compositional transformations, to coordinates
#clr
clr <- function(x, p=NULL){
  x <- format.check(x)
  log(x/geomeans.r(df = x, p = p))
}

#ilr
ilr <- function(x, V=NULL, p=NULL){
  x <- format.check(x)
  if (is.null(p)){
    p <- rep(1, ncol(x))
    names(p) <- colnames(x)
  }
  
  if (is.null(V)){ 
    V <- ilr_basis(ncol(x))
  }
  p <- p[colnames(x)]
  clr(x, p) %*% diag(p) %*% V
}

#transformation basis generation
#ilr base
ilr_basis <- function(n){
  ch <- stats::contr.helmert(n)
  ortho <- function(x) x/sqrt(sum(x^2))
  apply(ch, 2, FUN = ortho)
}

#transforms signary matrix to ILR basis, vector-wise
#when used with apply() outputs an ILR basis matrix with the same dimensions as input matrix
single.balance.basis <- function(x, p){
  pos <- names(x[which( x > 0)])
  neg <- names(x[which( x < 0)])
  if (isTRUE(length(pos)==0) | isTRUE(length(neg)==0)){
    x <- rep(NaN, length(x))
  } else {
    x <- replace(x, pos, 1/(length(pos)*p[pos]))
    x <- replace(x, neg, -1/(length(neg)*p[neg]))
    x <- x/sqrt(sum(x^2))
    x
  }
}

#application of single.balance.basis column-wise across input matrix. 
balance.basis <- function(W, p){
  apply(W, 2, single.balance.basis, p=p)
}

#reference shifts, as in Egozcue 2016
shift.reference <- function(x, p){
  x <- format.check(x)
  #Egozcue defines a shift in the reference measure on the simplex as a simple perturbation difference.
  df.shifted <- x/t(replicate(nrow(x), p))
  df.shifted
}

#create a predefined compositional balance with specific input taxa
create.balance <- function(df, num.tax, den.tax, weighted=FALSE, p=NULL){
  df <- format.check(df)
  
  #weights
  # throw error if unweighted and given p
  if (weighted==FALSE &&
      !is.null(p)){
    stop('If providing weights, please set "weighted=TRUE"')
  }
  
  if (weighted==FALSE){
    message('Using uniform weights (p)')
    p <- rep(1, ncol(df))
    names(p) <- colnames(df)
  }
  
  if (weighted==TRUE &&
      is.null(p)){
    message('Calculating part weights (p)')
    p <- geomeans.c(df)
  } else {
    p=p
  }
  p <- p[colnames(df)]
  
  #shift the reference measure by p
  message('Shifting the reference measure by p')
  df <- shift.reference(df, p)
  
  #create empty vector to transform to signary
  sigvec <- rep(0, ncol(df))
  names(sigvec) <- colnames(df)
  sigvec <- replace(sigvec, num.tax, 1)
  sigvec <- replace(sigvec, den.tax, -1)
  sigvec <- as.data.frame(sigvec)
  #perform ILR on generated balance
  V <- apply(sigvec, 2, single.balance.basis, p=p)
  ilr(x = df, V = V, p = p)
  
}
