## Auxiliary routines to compute all data for plotting the figures

## compute.fn must be a function that takes n.list, alpha, mu, eps and
## computes a data frame with all shadowed parameters and the 1-bit
## extractor seeds.
gen.plot.data <- function(n.list, alpha.list, mu.list, eps.list,
                         compute.fn) {
  data.n.mu.alpha.eps <- do.call(rbind, lapply(eps.list, function(.eps) {
    cat("eps: ", .eps, "\n")
    do.call(rbind, lapply(alpha.list, function(.alpha) {
      do.call(rbind, lapply(mu.list,
                            function(.mu) {
                              cat("  -> mu: ", .mu, "\n")
                              compute.fn(n.list, .alpha, .mu, .eps) })) })) }))
  
  data.n.mu.alpha.eps$m <- data.n.mu.alpha.eps$alpha*data.n.mu.alpha.eps$mu*data.n.mu.alpha.eps$n
  data.n.mu.alpha.eps$alpha <- as.factor(data.n.mu.alpha.eps$alpha)
  data.n.mu.alpha.eps$total.seed.2e <- data.n.mu.alpha.eps$seed.length**2
  data.n.mu.alpha.eps$seed.ratio <- data.n.mu.alpha.eps$n/(data.n.mu.alpha.eps$total.seed.2e)
  data.n.mu.alpha.eps$seed.inv.ratio <- 1/data.n.mu.alpha.eps$seed.ratio
  data.n.mu.alpha.eps$seed.ratio.m <- data.n.mu.alpha.eps$m/(data.n.mu.alpha.eps$total.seed.2e)
  data.n.mu.alpha.eps$seed.inv.ratio.m <- 1/data.n.mu.alpha.eps$seed.ratio.m
  data.n.mu.alpha.eps$extracted.bits <- data.n.mu.alpha.eps$m - data.n.mu.alpha.eps$total.seed.2e

  data.n.mu.alpha.eps$mu <- as.factor(data.n.mu.alpha.eps$mu)
  data.n.mu.alpha.eps$eps <- as.factor(data.n.mu.alpha.eps$eps)

  return(data.n.mu.alpha.eps)
}



## Compute the break-even point from which onwards we produce more randomness
## than we require for the seed.
## Input: data matrix produced by gen.plot.data where alpha does _NOT_ vary.
compute.breaks.alpha <- function(data.n.mu.eps, eps.list)  {
  dat.breaks <- do.call(rbind, lapply(eps.list, function(.eps) {
    data.frame(mu=mu.list, eps=.eps,
               break.even=sapply(mu.list, 
                 function(.mu) {
                   cat("eps: " , .eps, ", mu: ", .mu, "\n")
                   dat.subs <- subset(data.n.mu.eps, mu==.mu)
                   dat.subs <- subset(dat.subs, eps==.eps)

                   idx <- sum(dat.subs$seed.inv.ratio.m > 1)+1
                   return(dat.subs[idx,]$n)
                 })) }))

  return(dat.breaks)
}
