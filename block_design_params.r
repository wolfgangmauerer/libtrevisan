# Compute parameters for the combined weak design
# WM 09-Aug-2012
source("parameters.r")
source("blockdes.r")
library(ggplot2)

########### Check parameter behaviour ########### 
compute.block.params <- function(.m, .r, .t.seq) {
  dat <- data.frame(t=t.seq, m=.m, r=.r, a=sapply(.t.seq,
                                           function(t) { l.block(.r, .m, t); }))
  ## Total seed length is d=a*t^2 for the block design
  dat$total.seed <- sapply(dat$t, function(.t) { sub <- subset(dat, t == .t);
                                                 return(sub$a*(sub$t)**2) })
  ## Total seed length is d^2 for the weak design
  dat$total.seed.basic <- dat$t**2

  return(dat)
}

m.exp <- seq(4,8)
m.list <- sapply(m.exp, function(x) 10**x)
t.seq <- seq(8, 1024, 1)

dat.t <- do.call(rbind, lapply(m.list,
                               function(.m) compute.block.params(.m, 2*exp(1), t.seq)))
dat.t$overhead <- dat.t$total.seed/dat.t$total.seed.basic
dat.t$m <- as.factor(dat.t$m)

dat.tot <- data.frame(type="Block weak design", m=dat.t$m, total.seed=dat.t$total.seed)
dat.tot <- rbind(dat.tot,
                 data.frame(type="Standard weak design", m=dat.t$m, total.seed=dat.t$total.seed.basic))

g.sub <- ggplot(data=dat.t, aes(x=t, y=a, colour=m)) + geom_line() +
  scale_x_continuous("", breaks=(c(8,128,256,512,1024))) +
  scale_colour_discrete(name="Output length m") +
  ylab("Number of blocks") +
  opts(legend.position="none", plot.background=theme_rect(fill="white"))

g <- ggplot(data=dat.t, aes(x=t, y=total.seed, colour=m)) + geom_line() +
  scale_x_continuous("1-Bit extractor seed length",
                     breaks=(c(8,64, 128,256,512,1024))) +
  scale_y_log10("Block weak design seed length (bits) [log. scale]") +
  scale_colour_discrete(name="m (extracted\n # of bits)",
                        labels=math_format(10^.x)(m.exp),
                        breaks=m.list)

vp <- viewport(x=0.36, y=0.11, w=0.45, h=0.45, just=c("left", "bottom"))
pdf("paper/pictures/block_design_overview.pdf")
print(g)
print(g.sub, vp=vp)
dev.off()



###### Not relevant for plotting #######
dat.mi.leq.t <- data.frame(t=t.seq, m.leq.t=sapply(t.seq,
                                     function(t) { sum(m.list.block(2*exp(1), 100000, t) <= t); }))

## Distribution of the m_i
length(
       m.list.block(2*exp(1), 571787, 64)
       )
as.initialiser(2*exp(1), 571787, 64)

if (sum(m.list.block(2*exp(1), 100000, 2500)) != 100000) {
  cat("Hmmm..... sum(m_i) != m")
}
