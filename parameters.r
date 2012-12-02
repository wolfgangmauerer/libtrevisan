### Helper routines to determine parameters for 1-bit extractors

### Generic stuff
h <- function(x) {
  if (x < 2*.Machine$double.eps) {
    return(0);
  }
  if (x > 1-2*.Machine$double.eps) {
    return(0);
  }
  
  return(-x*log2(x) -(1-x)*log2(1-x))
}

delta <- function(eps, m) {
  return((eps/(3*sqrt(2)*m))**2)
}

### Auxiliary routines for the XOR extractor ###

## We minimise |gamma - (l, delta)| to obtain the proper value 
entrop.argument <- function(l, eps) log(2)/l*log2(2/eps)

## Lower bound on l
l.lower.xor <- function(eps) {
  return(log(2)/0.5*log(2/eps))
}

## Note: The following definitions hold
## k = alpha*n (i.e., the min-entropy produced by the source)
## m = mu*k    (number of extracted bits. mu is the extraction fraction, usually
##              small -- say, 1/10 or even less)

gamma.func <- function(alpha,mu,r,n,eps) {
  return(((alpha*n)-r*mu*alpha*n - 6*log2((1+sqrt(2))/eps) - log2(4/3))/n)
}

opt.func.xor <- function(l, gamma, eps) {
  if (entrop.argument(l, eps) > 1 || entrop.argument(l, eps) < 0) {
    stop("Parameter out of range -- binary entropy cannot be computed for p < 0 or p > 1")
  }

  return(abs(gamma - h(entrop.argument(l, eps))))
}

## do.opt.xor delivers the optimal l for the given parameters
do.opt.xor <- function(alpha,mu,r,n,eps) {
  l.opt <- optimize(function(x) opt.func.xor(x, gamma.func(alpha,mu,r,n,eps), eps),
                    interval=c(l.lower.xor(eps), 8000))$minimum
  l.opt <- ceiling(l.opt)
  return(l.opt)
}

## ########### Auxiliary routines for the Reed-Solomon-Hadamard code #############
## t = 2*l
l.rsh <- function(n, eps) {
  return (ceiling(log2(n)+2*log(2/eps)))
}

r.rsh <- function(n, eps) {
  return (ceiling(n/l.rsh(n, eps)))
}



## ########### Auxiliary routines for the Lu extractor code #############
lu.nu <- function(lambda, delta, l) {
  return(1+lambda**2-delta**(4/l))
}

## We want to solve w log(w) =(1-nu+w)log(1-nu+w), and can do this
## by computing the minimum of |wlog(w) - () log()|
opt.func.lu <- function(w, nu) {
  return(abs(w*log2(w)-(1-nu+w)*log2(1-nu+w)))
}

c.lu <- function(w, lambda0) {
  return(ceiling(log2(w)/(2*log2(lambda0))))
}

l.lu <- function(eps, m, nu, w) {
  return(ceiling(4*log2(delta(eps,m))/log(1-nu+w)))
}

## Given nu and lambda (implicitly fixed by the graph
## construction), compute the optimal value of w, and derive c from it
do.opt.lu <- function(nu) {
  w <- optimize(function(x) opt.func.lu(x, nu), interval=c(1e-5, 1))$minimum

  return (w)
}

## Given a set of initial values (experimental constraints, and a chosen
## nu), determine w, and infer c and l)
do.compute.lu <- function(nu, m, eps, lambda0) {
    w <- do.opt.lu(nu)
    c <- c.lu(w, lambda0)
    l <- l.lu(eps, m, nu, w)

    return(data.frame(w=w, c=c, l=l))
}

steps.lu <- function(nu, m, eps, lambda0) {
  res <- do.compute.lu(nu, m, eps, lambda0)
  return (res$c*res$l)
}

k.lu <- function(.nu, .n, .m, .r, .eps) {
  return(h(.nu)*.n + .r*.m + 6.0*log((2.0+sqrt(2.0))/.eps) - 2.0)
}

seed.lu <- function(n, c, l) {
    return(log2(n) + 3*c*(l-1) + l)
}
