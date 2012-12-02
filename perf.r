### Plot performance measurement data
library(ggplot2)
library(scales)

plot.perf.per.core <- function(.dat, .colour=NA) {
  p <- ggplot(data=.dat, aes(x=cpus, y=perf.per.core)) +
    scale_y_continuous("Throughput [kbit/s] per core") +
    xlab("Number of CPU Cores")

  if (!is.na(.colour)) {
    p <- p + geom_boxplot(aes(x=factor(cpus), colour=.colour))
  } else {
    p <- p + geom_boxplot(aes(x=factor(cpus)))
  }
  
  return(p)
}

plot.perf <- function(.dat, .colour=NA) {
  p <- ggplot(data=.dat, aes(x=cpus, y=perf)) +
    scale_y_continuous("Throughput [kbit/s]") +
    xlab("Number of CPU Cores")

  if (!is.na(.colour)) {
    p <- p + geom_boxplot(aes(x=factor(cpus), colour=.colour))
  } else {
    p <- p + geom_boxplot(aes(x=factor(cpus)))
  }
  
  return(p)
}

plot.rsh.block.scale.n <- function(.dat) {
  p <- ggplot(data=.dat, aes(x=n, y=perf, colour=wd)) +
    geom_boxplot(aes(x=as.factor(n))) +
    scale_y_continuous("Throughput [kbit/s]") +
    scale_x_discrete("Number of input bits (log. scale)",
                     labels=math_format(2^.x)(dat.block.rsh.mb$n),
                     breaks=dat.block.rsh.mb$n) +
                       scale_colour_discrete(name="Weak design")
  
  return(p)
}



######################
# Determine how block(gf{2x,p}) and the rsh extractor scale with varying input length
# (measurement for per-cpu and per-n scaling is in one file the the MacBook)
dat.block.rsh.mb <-
  read.table("meas/perf_block_rsh_macbook_n.txt", header=T) # RSH measurement on MacBook Air
dat.block.rsh.mb$perf.per.core <- dat.block.rsh.mb$perf/dat.block.rsh.mb$cpus

dat.block.rsh.opteron <-
  read.table("meas/perf_block_rsh_opteron_n.txt", header=T) # RSH measurement on 48 opteron
dat.block.rsh.opteron$perf.per.core <- dat.block.rsh.opteron$perf/dat.block.rsh.opteron$cpus

dat.ma <- read.table("meas/perf_ma.txt", header=T)
dat.ma$cpus=1
dat.ma$perf <- dat.ma$perf/1000 # Ma et al. specify their results in bits/s

# TODO: Measure with more iterations
#p <- plot.rsh.block.scale.n(subset(dat.block.rsh.mb, cpus==4))
#p <- p + opts(title="RSH extractor throughput (notebook, 4 cores, 10 repetitions)")
##print(p)
#ggsave("paper/pictures/rsh_designs_macbook_scale_n.pdf", p)

#p <- plot.rsh.block.scale.n(dat.block.rsh.opteron)
#p <- p + opts(title="RSH extractor throughput (workstation, 48 cores, 30 repetitions)")
##print(p)
#ggsave("paper/pictures/rsh_designs_opteron.pdf", p)

dat.block.rsh.mb$type <- "MacBook (4 cores, 40 repetitions)"
dat.block.rsh.opteron$type <- "Opteron (48 cores, 30 repetitions)"

p <- plot.rsh.block.scale.n(rbind(dat.block.rsh.mb, dat.block.rsh.opteron)) +
  facet_grid(type~., scale="free_x")
#print(p)
ggsave("paper/pictures/rsh_designs_scale_n.pdf", p)


## Compare with Ma et al.'s data
p <- ggplot(data=dat.block.rsh.mb[dat.block.rsh.mb$cpus %in% c(1,4) &
              dat.block.rsh.mb$wd=="block(gf2m)",],
            aes(x=as.factor(n), y=perf)) +
  geom_boxplot(aes(fill=as.factor(cpus))) + 
  scale_x_discrete("Number of input bits",
                   labels=math_format(2^.x)(dat.block.rsh.mb$n),
                   breaks=dat.block.rsh.mb$n) +
  geom_point(data=dat.ma, shape=17, size=2.5) +
  scale_y_log10("Throughput [kbit/s] (log. scale)") +
  scale_fill_discrete("Number of cores") +
  opts(title=expression(paste("Throughput scaling (block(GF(", 2^x, "), RSH) (40 iterations)")))
#print(p)
ggsave("paper/pictures/rsh_ma_compare.pdf", p)

## Scaling behaviour of same primitives for fixed n and varying number of CPUs
## TODO: Use this figure as inset in the opteron graph below.
## TODO: Note that the MacBook outperforms the Opteron by a factor of two for
## the single-CPU power.
dat.block.rsh.mb$cpus <- ordered(dat.block.rsh.mb$cpus, levels=1:4)
p <- ggplot(data=subset(dat.block.rsh.mb, n==16), aes(x=cpus, y=perf.per.core, colour=wd)) +
  geom_boxplot() +
  scale_y_continuous("Throughput [kbit/s]") + xlab("Number of CPU Cores")
print(p)

##################
## Scaling behaviour of block(gfp)+rsh on the 48 core opteron
dat.rsh.opteron <- read.table("meas/perf_block_rsh_opteron_cpus.txt", header=T)
## Plot with factors 1,2,4,6,... will have incorrect distances otherwise
# TODO: Use an ordered factor, and introduce empty elements for the unused elements.
# This restores te distance scaling
#dat.rsh.opteron <- subset(dat.rsh.opteron, cpus > 1)
dat.rsh.opteron$i <- as.factor(dat.rsh.opteron$i)
dat.rsh.opteron$n <- as.factor(dat.rsh.opteron$n)
dat.rsh.opteron$perf.per.core <- dat.rsh.opteron$perf/dat.rsh.opteron$cpus

## Fit a mixed effects model to check the speedup behaviour
#fm <- lmer(perf~cpus+(1|i), dat.rsh.opteron)
#dat.rsh.opteron$perf.mlm <- coef(fm)$i[1,1]+dat.rsh.opteron$cpus*coef(fm)$i[1,2]

# Does not really work because of the factor scaling...
#p <- ggplot(data=dat.rsh.opteron, aes(x=as.factor(cpus), y=perf)) + geom_boxplot() +
#  scale_y_continuous("Throughput [kbit/s]") + xlab("Number of CPU Cores") +
#  geom_point(aes(y=perf.mlm), colour="red")
#print(p)

p <- plot.perf(dat.rsh.opteron) +
  opts(title=expression(paste("Scaling behaviour: Block(GF(p))+RSH on 48 core Opteron (n=",
      2**16, ", 40 repetitions)")))
#print(p)
ggsave("paper/pictures/block_rsh_cores_scaling.pdf", p)

## XOR+gfp performance for varying input sizes n (on Opteron with 48 cores)
## m is fixed to n/100
dat.xor <- read.table("meas/xor_gfp_opteron_scale_n.txt", header=T)
dat.xor$n <- dat.xor$n/1000 # Convert to GiBit (dividing by 1000 instead of 1024
                            # gives nicer numbers)
p <- ggplot(data=dat.xor, aes(x=n, y=perf)) + geom_boxplot(aes(x=factor(n))) +
    scale_y_continuous("Throughput [kbit/s]") +
    xlab("Input bit length n [GiBit], m=n/100") +
  opts(title="XOR+GF(p) throughput scaling (48 cores, 40 repetitions)")
#print(p)
ggsave("paper/pictures/xor_gfp_n_scaling.pdf", p)


## XOR+gfp performance scaling with varying number of cores (on opteron)
dat.xor.cores <- read.table("meas/xor_scaling_opteron_cpus.txt", header=T)
dat.xor.cores$perf.per.core <- dat.xor.cores$perf/dat.xor.cores$cpus

#p <- plot.perf(dat.xor.cores)
#print(p + xlab("Number of CPU cores (n=200MiBit, m=n/100)"))

#p <- plot.perf.per.core(dat.xor.cores)
#print(p + xlab("Number of CPU cores (n=200MiBit, m=n/100)"))

## Combined figure for xor and rsh for per-core scaling
## xor: n=200MiBit, m=n/100
## rsh: n=2**16, m=2**15
dat.comb <- rbind(cbind(dat.xor.cores, type="XOR+GF(p) (40 iterations)"),
                  cbind(dat.rsh.opteron, type="RSH+Block(GF(p)) (40 iterations)"))
p <- plot.perf.per.core(dat.comb) + facet_grid(type~., scales="free_y") +
  opts(title=expression(paste("XOR: n=", 2%*%10^6, " (m=n/100); RSH: n=", 2^16,
      ", m=", 2^15, " (48 cores)")))
#print(p)
ggsave("paper/pictures/xor_rsh_per_core_comparison.pdf", p)
