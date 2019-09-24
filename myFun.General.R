get.local.lib <- function(root = './'){
  local.lib = paste(root, 'Rlib', sep = '/')
  .libPaths(local.lib)
  return(local.lib)
}

is.installed <- function(pkg, local.lib = get.local.lib()){
  return(pkg %in% rownames(installed.packages(lib.loc = local.lib)))
}

mat.fill.row <- function(mat, row.names = rownames(mat), default=0){
  is.vector = F
  if (!is.character(row.names)){
    row.names = as.character(row.names)
  }
  if (is.vector(mat)){ # 20140629
    mat = matrix(mat, nrow=length(mat), ncol=1, dimnames=list(names(mat), 'mat'))
    is.vector = T
  }
  mat.new = matrix(default, nrow=length(row.names), ncol=ncol(mat), dimnames=list(row.names, colnames(mat)))
  i = match(row.names, rownames(mat))
  i.i = !is.na(i)
  mat.new[i.i,] = mat[i[i.i],]
  if (is.vector)
    mat.new = mat.new[,1]
  return(mat.new)
}


cbind.union <- function(x, y, default=0){
  # union the rownames and combine columns
  # Yong Fuga Li
  
  if (is.vector(x)){
    x = matrix(x, nrow=length(x), ncol=1, dimnames=list(names(x), 'x'))
  }
  if (is.vector(y)){
    y = matrix(y, nrow=length(y), ncol=1, dimnames=list(names(y), 'y'))
  }
  if (is.null(x))
    return(y)
  if (is.null(y))
    return(x)  
  names.1 = rownames(x)
  if (is.null(names.1))
    stop('Input x has no names or rownames')
  
  names.2 = rownames(y)
  if (is.null(names.2))
    stop('Input y has no names or rownames')
  
  names.all = union(names.1, names.2);
  x = mat.fill.row(x, names.all, default = default)
  y = mat.fill.row(y, names.all, default = default)
  return(cbind(x,y))
}


max2 <- function(x, y){
  if (length(x)>1 & length(y)==1){
    y = x*0 + y
  }else if (length(x)==1 & length(y)>1){
    x = y * 0 + x
  }else if (length(x)!=length(y)){
    stop('x and y are of unequal sizes')
  }
  i = x<y
  x[i] = y[i]
  return(x)
}

min2 <- function(x, y){
  return(- max2(- x, - y))
}

seq2 <- function(...){
  tryCatch(
    return(seq(...)),
    error=function(e){
      return(c())
    },
    finally={
    }
  )
}

changenames <- function(fromRe = '^nack2Mix5_B3DTF(.*)$', toRe ='nack2Mix2_B3DTF\\1'){
  files = list.files(pattern = fromRe)
  for (f in files){
    newf = sub(f, pattern = fromRe, replacement = toRe)
    system(paste('mv', f, newf))
  }
}


logsig <- function(x) 1/(1+exp(-x))
softmax <- function(x) { y = exp(c(x,0)); return(y/sum(y))} # generalie logsig
softmax.neg <- function(x) { y = exp(c(x,0)); return((sum(y)-y)/sum(y))} # generalie logsig

gradient.softmax <- function(x){
  p = softmax(x)
  n = length(p)
  q = vector(mode = 'numeric', n)
  for (i in 1:n){
    q[i] = sum(p[setdiff(1:n, i)])
  }
  n = length(p)
  grad = diag(p, nrow = n) - p %*% t(p)
  diag(grad) = p * q # to avoid numerical issues of p close to 1 
  grad[1:n, 1:(n-1)] # note the last x is fixed as 0, we do not the gradient for it.
}

logit <- function(x) return(log(x/(1-x)))

regexpr.match <- function(pat, txt, perl=T,...){
  all.matches = c()
  for (t in txt){
    aa <- regexpr(pat,t, perl=perl,...)
    st <- attr(aa,'capture.start')
    en <- attr(aa,'capture.start')+attr(aa,'capture.length')-1
    matches <- c()
    for (i in 1:length(st)){
      matches[i] <- substr(t,st[i],en[i])      
    }
    all.matches = rbind(all.matches, matches)
  }
  
  return(all.matches)    
}

rowLogSumExp <- function(x){
  # 20160919, YFL
  m = rowMax(x)
  return(log2(rowSums(2^(x-m)))+m)
}

rowMax <- function(x,...){
  if (all(is.na(x))){
    return(NA)
  }else{
    return(apply(x,1,max,...))    
  }
}