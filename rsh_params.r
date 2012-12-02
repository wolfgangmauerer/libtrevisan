## Parameter calculations for the Reed-Solomon-Hadamard code
source("parameters.r")
source("plot_calculations.r")
library(ggplot2)
library(scales)

n.list <- sapply(seq(4, 10, length.out=1500), function(x) 10**x)
## RSH does not depend on alpha and eps
alpha.list <- c(1)
mu.list <- c(1)
eps.list <- c(10**-3, 10**-7, 10**-11, 10**-15)

compute.n.rsh <- function(n.list, .alpha, .mu, .eps) {
  return(data.frame(n=n.list, alpha=.alpha, mu =.mu, eps=.eps,
                     seed.length=sapply(n.list, function(n) 2*l.rsh(n, .eps)),
                     poly.degree = sapply(n.list, function(n) r.rsh(n, .eps))))
}

dat <- gen.plot.data(n.list, alpha.list, mu.list, eps.list, compute.n.rsh)

g.sub <- ggplot(data=dat, aes(x=n, y=seed.length, colour=eps)) +
  geom_line() +
  scale_x_log10("", breaks=c(10**4, 10**6, 10**8, 10**10, 10**12)) +
  scale_y_continuous("Seed t") +
  opts(legend.position="none", plot.background=theme_rect(fill="white"))

g.sub2 <- ggplot(data=dat, aes(x=n, y=seed.inv.ratio.m, colour=eps)) +
  geom_hline(aes(yintercept=1), colour="red", size=1, linetype="dashed") + geom_line() +
  scale_x_log10("", breaks=c(10**4, 10**6, 10**8, 10**10, 10**12)) +
  scale_y_log10("Seed length/m [log. scale]") +
  opts(legend.position="none", plot.background=theme_rect(fill="white"))

g <- ggplot(data=dat, aes(x=n, y=poly.degree, colour=eps)) + geom_line() +
  scale_x_log10("n (# of input bits) [log. scale]",
                breaks=c(10**4, 10**6, 10**8, 10**10, 10**12)) +
  scale_y_log10("Polynomial degree [log. scale]") +
  scale_colour_discrete(name=expression(paste("Bit error ", epsilon)),
                        labels=math_format(10^.x)(log10(eps.list)),
                        breaks=eps.list) +
  opts(title="RSH extractor: Scaling behaviour")

vp <- viewport(x=0.82, y=0.1, w=0.35, h=0.35, just=c("right", "bottom"))
vp2 <- viewport(x=0.12, y=0.56, w=0.35, h=0.35, just=c("left", "bottom"))
pdf(file="paper/pictures/rsh_overview.pdf")
print(g)
print(g.sub, vp=vp)
print(g.sub2, vp=vp2)
dev.off()
