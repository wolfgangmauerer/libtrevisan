## ######### Auxiliary routines for the block weak design ############
## r is the r overlap value for the basic design
## m denotes the output length
## t is the seed length required by the 1-bit extractor
l.block <- function(r, m, t) {
  res <- ceiling(log2((m-r)/(t-r))/log2(r/(r-1)));
  if (t > m)
    stop("Bogus parameters: t > m")
  if (r <= 0 || m <= 0)
    stop("Bogus parameter: r <= 0 or m <= 0")
  
  return(max(1, res));
}

## i is the number of the sub-design, 0 <= i < l
## r is the r overlap value for the basic design
## m denotes the output length
n_i.block <- function(i, r, m) {
  return(((1-1/r)**i)*(m/r-1))
}

## Compute a list of {m_i}, that is, the number of lines in the ith
## design. Parameters as above. A recursive implementation turned out
## to be very slow, so changed to a more boring side-effect scheme.
m.list.block <- function(r, m, t) {
  l <- l.block(r,m,t)
  m.seq <- rep(0, l)
  n.seq <- rep(0, l)
  
  for (j in 0:(l-1)) {
    n.seq[j+1] <- n_i.block(j, r, m)
    m.seq[j+1] <- ceiling(sum(n.seq[1:(j+1)]))
    if (j > 0) {
      m.seq[j+1] <- m.seq[j+1] - sum(m.seq[1:j])
    }
  }

  m.seq[l+1] <- m - sum(m.seq)
  return(m.seq)
}

## Take the information from m.list.block to produce the block descriptor syntax in form
## of a C++ array initialiser
## (The routine assumes that lst[1] > lst[length(lst)])
compute.descriptor <- function(lst, desc, consumed, fixup) {
  if(length(lst) == 1) {
    return(rbind(desc, data.frame(highest=consumed+lst[1], iter=lst[1],
                                  bits.per.iter=1, fixup=F)))
  }

  next.fixup = FALSE
  if(lst[length(lst)-1] < lst[length(lst)]) {
    next.fixup = TRUE
    iter <- lst[length(lst)-1]
  } else {
    iter <- lst[length(lst)]
  }

  bits.per.iter <- length(lst)
  consumed = consumed + iter*bits.per.iter

  if (is.null(desc)) {
    desc <- data.frame(highest=consumed, iter=iter,
                       bits.per.iter=bits.per.iter, fixup=fixup)
  } else {
    desc <- rbind(desc, data.frame(highest=consumed, iter=iter,
                                   bits.per.iter=bits.per.iter, fixup=fixup))
  }

  lst  <- lst - iter
  compute.descriptor(lst[lst > 0], desc, consumed, next.fixup)
}

as.initialiser <- function(r, m, t) {
  lst <- m.list.block(r, m, t)

  ## Get the recursion in compute.descriptor going
  df <- compute.descriptor(lst, NULL, 0, FALSE)

  p.i <- function(i) {
    return(sprintf("%d", i))
  }

  to.cpp <- function(x) {
      return(as.character(paste(c("{ ", p.i(x[1]), ",", p.i(x[2]), ",",
                                  p.i(x[3]), ",", p.i(x[4]) , " } "),
                                collapse="")))
    }

  ## Turn the data frame into a C++ array initialiser
  res <- paste("{", paste(apply(df, 1, to.cpp), collapse=", ", sep=""), "}", collapse=" ")

  return(res)
}

block.param.as.matrix <- function(r, m, t) {
  lst <- m.list.block(r, m, t)

  ## Get the recursion in compute.descriptor going
  df <- compute.descriptor(lst, NULL, 0, FALSE)

  return(as.matrix(df))
}
