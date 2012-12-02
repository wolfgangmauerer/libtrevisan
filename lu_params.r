## Optimisation script for the Lu extractor
## (C) 2012, Wolfgang Mauerer <wm@linux-kernel.net>
## Licensed under the GPLv2

source("parameters.r")
source("plot_calculations.r")
library(ggplot2)
library(scales)

## Result of do.compute.lu: data.frame with w, c, l
## Number of random walk steps: c*l
#eps <- 1e-7
#m <- 2**20
#nu <- 0.49
## do.compute.lu(nu, m, eps, lambda0)
#test <- do.compute.lu(0.45, 2**15, eps, lambda0)
#cat("Computed c=", test[[2]], ", l=", test[[3]], ", w=", test[[1]],
#    ", nu=", nu, "\n", sep="")
#
#cat("Derived value of nu: ", lu.nu(lambda0**test[[2]], delta(eps, m), l), "\n", sep="")
#cat("Derived value of lambda^2: ", (lambda0**test[[2]])**2, "\n\n", sep="")

nu.list <- c(0.1, 0.2, 0.3, 0.4, 0.45)
eps.base <- c(-3, -7, -15)
eps.list <- sapply(eps.base, function(x) 10**x)
m.list <- seq(10**4, 10**9, length.out=100)
lambda0 <- 5*sqrt(2)/8

dat.lu.steps <- do.call(rbind, lapply(nu.list, function(.nu) {
  do.call(rbind, lapply(eps.list, function(.eps) {
    do.call(rbind, lapply(m.list, function(.m) {
      return(data.frame(nu=.nu, m=.m, eps=.eps, steps=steps.lu(.nu, .m, .eps, lambda0)))
    }))
  }))
}))

g <- ggplot(dat=dat.lu.steps, aes(x=m, y=steps)) + geom_line(aes(colour=as.factor(nu))) +
  facet_grid(~eps) + scale_y_log10("# of random walk steps [log. scale]") +
  scale_x_log10("# output bits (m) [log. scale]") +
  scale_colour_discrete(name=expression(nu)) +
  opts(title="Lu extractor: Number of random walk steps")
#print(g)
ggsave("paper/pictures/lu_steps.pdf", g)
# -> the number of steps grows considerably with decreasing mu. Since all except tiny
# extraction ratios require large values of nu, the extractor is effectively only
# usable for very small output lengths


## Determine the maximal value of mu (which leads to the smallest number
## of random walk steps) that allows for extracting mu*alpha*n bits from n bits
max.nu.lu <- function(n, mu, alpha) {
  ## We could also use a canned constrained nonlinear optimisation routine
  ## here, but the optimisation is simple enough to be performed by brute force
  ## (note: m = mu*alpha*n)
  nu.list <- seq(0.005, 0.49, length.out=100)
  res <-  do.call(rbind, lapply(nu.list, function(.nu) {
    data.frame(nu=.nu, val=alpha*n - k.lu(.nu, n, mu*alpha*n, 2*exp(1), 10**(-7)))
  }))

  return(max(res[res$val > 0,]$nu))
}

# NOTE: Calculation noes not depend on n
mu.list <- seq(10**(-5), 10**(-1), length.out=250)
alpha.list=c(0.8, 0.85, 0.9, 0.95, 0.99)
lu.steps.dat <- do.call(rbind, lapply(alpha.list, function(.alpha) {
  do.call(rbind,
          lapply(mu.list,
                 function(.mu) {
                   return(data.frame(mu=.mu, alpha=.alpha,
                                     nu=max.nu.lu(10**9, .mu, .alpha)))}
                 ))
}))

g <- ggplot(dat=lu.steps.dat, aes(x=mu, y=nu, colour=as.factor(alpha))) + geom_line() +
  scale_x_log10(expression(paste("Extraction factor ", mu))) +
  scale_y_continuous(expression(paste("Parameter ", nu))) +
  scale_colour_discrete(name=expression(paste("Entropy factor ", alpha))) +
  opts(title=expression(paste("Lu extractor: Interplay between ", nu, " and ", mu)))
#print(g)
ggsave("paper/pictures/lu_param_nu.pdf", g)

## Number of extracted bits vs. seed length plots (same as for the XOR extractor)
n.list <- sapply(seq(5, 12, length.out=30), function(x) 10**x)
alpha.list <- seq(0.7, 0.95, 0.05)
mu.list <- c(0.001, 0.05, 0.01, 0.1)
## Use the same parameter dispatcher as for the XOR extractor so that
## we can re-use the infrastructure. 
compute.n.lu <- function(n.list, .alpha, .mu, .eps) {
  return(data.frame(n=n.list, alpha=.alpha, mu=.mu, eps=.eps,
                    seed.length=sapply(n.list, function(n) {
                      nu <- max.nu.lu(n, .mu, .alpha)
                      m <- .alpha*.mu*n
                      lu.res <- do.compute.lu(nu, m, .eps, lambda0)
                      return(seed.lu(n, lu.res$c, lu.res$l))
                    })))
}


## Scaling behaviour. Mostly c&p from the XOR extractor plots.
data.n.mu.alpha.eps <- gen.plot.data(n.list, alpha.list, mu.list, eps.list, compute.n.lu)

g <- ggplot(data=data.n.mu.alpha.eps, aes(x=n, y=seed.length, colour=alpha)) + geom_line() +
  scale_x_log10("n (# of input bits) [log. scale]", breaks=c(10**5, 10**7, 10**9, 10**12)) +
  scale_y_continuous("Total seed length for r=2e weak design") +
  facet_grid(eps~mu, scales="free_y") +
  scale_colour_discrete(name=expression(paste("Source\nentropy ", alpha))) +
  opts(title="Lu extractor: Scaling behaviour")
print(g)
ggsave("paper/pictures/lu_overview.pdf", g)

dat.subs <- subset(data.n.mu.alpha.eps, eps==10**(-7))
dat.subs <- subset(dat.subs, mu==0.05)
g <- ggplot(data=dat.subs, aes(x=n, y=seed.inv.ratio.m, colour=alpha)) +
  geom_hline(aes(yintercept=1), colour="red", size=1, linetype="dashed") + geom_line() +
  scale_x_log10("n (# of input bits) [log. scale]") +
  scale_y_log10(expression(paste("Ratio of seed length to extracted bits (",
      m=mu*alpha*n, ") [log. scale]"))) +
  scale_colour_discrete(name=expression(paste("Source\nentropy ", alpha))) +
  opts(title=expression(paste("Lu extractor: Seed to extraction ratio (", mu, "=0.05)")))
#print(g)
ggsave("paper/pictures/lu_ratio.pdf", g)
