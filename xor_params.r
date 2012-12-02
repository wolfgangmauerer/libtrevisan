# Simple optimisation script for the XOR extractor
# (C) 2012, Wolfgang Mauerer <wm@linux-kernel.net>
# Licensed under the GPLv2
source("parameters.r")
source("plot_calculations.r")
library(ggplot2)
library(scales)

## Note: 10^9 is one Gigabit
## function(alpha,mu,r,n,eps) {
#do.opt.xor(0.9, 0.01, 2*exp(1), 10**9, 10**(-15))

seed.xor <- function(n) {
  l <- do.opt.xor(0.8, 0.05, 2*exp(1), n, 10**(-7))
  cat ("l: ", l)
  return(l*log2(n))
}

#seed.xor(10**10)
# Extractable bits (maximal m): 0.8*0.05 * 10^10 = 0.04 * 10^10 = 4*10^8


## Compute systematic parameter comparison graphs. The length of the
## XOR extractor seed is t=l*log2(n), which is what this loop computes
compute.n.xor <- function(n.list, .alpha, .mu, .eps) {
  return(data.frame(n=n.list, alpha=.alpha, mu =.mu, eps=.eps,
                    seed.length=sapply(n.list, function(n)
                      do.opt.xor(.alpha, .mu, 2*exp(1), n, .eps)*log2(n))))
}
#data.n <- compute.n.xor(n.list, 0.8, 0.06)
#ggplot(data=data.n, aes(x=n, y=seed.length)) + geom_line() + xlab("n (# of input bits)") +
#  ylab("1-bit extractor seed length")

n.list <- sapply(seq(5, 12, length.out=1500), function(x) 10**x)
alpha.list <- seq(0.7, 0.95, 0.05)
mu.list <- c(0.01, 0.1, 0.15, 0.17)
eps.list <- c(10**-3, 10**-7, 10**-15)

data.n.mu.alpha.eps <- gen.plot.data(n.list, alpha.list, mu.list, eps.list, compute.n.xor)

g <- ggplot(data=data.n.mu.alpha.eps, aes(x=n, y=seed.length, colour=alpha)) + geom_line() +
  scale_x_log10("n (# of input bits) [log. scale]", breaks=c(10**5, 10**8, 10**10, 10**12)) +
  scale_y_continuous("Seed length t for r=2e weak design") +
  facet_grid(eps~mu, scales="free_y") +
  scale_colour_discrete(name=expression(paste("Source\nentropy ", alpha))) +
  opts(title="XOR extractor: Scaling behaviour")
#print(g)
ggsave("paper/pictures/xor_overview.pdf", g)

## Not too informative because all panel elements have very similar content
g <- ggplot(data=data.n.mu.alpha.eps, aes(x=n, y=seed.inv.ratio.m, colour=alpha)) +
  geom_line() +
  scale_x_log10("n (# of input bits) [log. scale]", breaks=c(10**5, 10**8, 10**10, 10**12)) +
  scale_y_log10("Ratio of seed length to extracted bits [log. scale]") +
  facet_grid(eps~mu, scales="free_y")
#print(g)
#ggsave("paper/pictures/xor_ratio_panel.pdf", g)

dat.subs <- subset(data.n.mu.alpha.eps, eps==10**(-7))
dat.subs <- subset(dat.subs, mu==0.1)
g <- ggplot(data=dat.subs, aes(x=n, y=seed.inv.ratio.m, colour=alpha)) +
  geom_hline(aes(yintercept=1), colour="red", size=1, linetype="dashed") + geom_line() +
  scale_x_log10("n (# of input bits) [log. scale]") +
  scale_y_log10(expression(paste("Ratio of seed length to extracted bits (",
      m=mu*alpha*n, ") [log. scale]"))) +
  scale_colour_discrete(name=expression(paste("Source\nentropy ", alpha))) +
  opts(title="XOR extractor: Seed to extraction ratio")
#print(g)

dat.subs$extracted.bits[dat.subs$extracted.bits <= 0] <- NA 
g.sub <- ggplot(data=dat.subs, aes(x=n, y=extracted.bits, colour=alpha)) +
  geom_line() +
  scale_x_log10("", breaks=c(10**5, 10**8, 10**11)) +
  scale_y_log10("m-d") +
  facet_grid(~mu, scales="free_y") +
  opts(legend.position="none", plot.background=theme_rect(fill="white"))
vp <- viewport(x=0.82, y=0.56, w=0.35, h=0.35, just=c("right", "bottom"))
pdf(file="paper/pictures/xor_ratio.pdf")
print(g)
print(g.sub, vp=vp)
dev.off()


## Find the break-even points numerically for a range of mu and eps (for given alpha=0.8)
mu.list <- seq(0.005, 0.15, 0.005) # Note: We run into numerical problems for > 0.18
eps.base <- seq(-15, -7, 2)
eps.list <- sapply(eps.base, function(x) 10**x)
data.n.mu.eps <- gen.plot.data(n.list, c(0.8), mu.list, eps.list, compute.n.xor)
dat.breaks <- compute.breaks.alpha(data.n.mu.eps, eps.list)

g <- ggplot(data=dat.breaks, aes(x=mu, y=break.even, color=as.factor(eps))) +
  geom_line() +
  scale_x_continuous(expression(paste("Extraction fraction ", mu)),
                     breaks=c(0.005, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14)) +
  scale_y_log10("n at break-even point [log. scale]") +
  opts(title="XOR extractor: Break-even points") +
  scale_colour_discrete(name=expression(paste("Bit error ", epsilon)),
                        labels=math_format(10^.x)(eps.base),
                        breaks=eps.list)
print(g)
ggsave("paper/pictures/xor_break_even.pdf", g)


########## Old stuff that is not required for the paper #####
ggplot(data=data.n.alpha, aes(x=n, y=total.seed, colour=as.factor(alpha))) + geom_line() +
  xlab("n (# of input bits)") + ylab("Total seed length for r=2e weak design")

g <- ggplot(data=data.n.alpha, aes(x=n, y=seed.ratio, colour=as.factor(alpha))) +
  geom_line() + xlab("n (# of input bits)") + ylab("Ratio input length to total seed")
#print(g)
#ggsave("/tmp/xor.pdf", g)
