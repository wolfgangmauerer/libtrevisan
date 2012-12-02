library(ggplot2)
dat <- read.table("params/params_curr.txt", header=TRUE)
# Eliminate the cases where maxsteps was reached
dat <- dat[which(dat$mu <= 0.75),]

# The number of random walk steps is essentially the same irregardless
# the degree, but goes up extremely fast for larger mu
dat$deg <- as.factor(dat$deg)
plt <- ggplot(data=dat, aes(x=mu, y=zeta)) + geom_point(aes(colour=deg)) +
  geom_line(aes(colour=deg)) +
  opts(title="Random walk steps dependent on degree") +
       labs(y="zeta [log scale]", x="mu") + scale_y_log10()
#print(plt)
ggsave("paper/pictures/params_curr.pdf", plt)

# The optimal value for lambda brings decisive computational advantages
# The value for d, in contrast, is not so important
dat <- read.table("params/params_opt.txt", header=TRUE)
dat$deg <- as.factor(dat$deg)
plt <- ggplot(data=dat, aes(x=mu, y=zeta)) + geom_point(aes(color=deg)) +
  facet_wrap(~d, scales="free") + 
  opts(title="Random walk steps (facets: d) for optimal lambda") +
       labs(y="zeta", x="mu") #+ scale_y_log10()
#print(plt)
ggsave("paper/pictures/params_opt.pdf", plt)


# Fix d=8. How does lambda influence the results when we go
# travese the lambda range from optimal to currently possible?
dat <- read.table("params/params_d=8_lambda.txt", header=TRUE)
dat <- dat[which(dat$deg==2),]
dat <- dat[which(dat$mu <= 0.75),]
dat$lambda <- as.factor(dat$lambda)
plt <- ggplot(data=dat, aes(x=mu, y=zeta)) + geom_point(aes(colour=lambda)) +
  geom_line(aes(colour=lambda)) + scale_y_log10() +
  opts(title="Random walk steps for d=8 and varying lambda") +
       labs(y="zeta [log scale]", x="mu") #+ scale_y_log10()
ggsave("paper/pictures/extender_d8_lambda.pdf", plt)

# Compare the performance of weak design and 1-bit-extractor for different degrees
dat.wd <- data.frame(deg=2, read.table("params/stat_deg2/wd_stats.txt", col.names=c("val")))
dat.wd <- rbind(dat.wd,
                data.frame(deg=3, read.table("params/stat_deg3/wd_stats.txt",
                             col.names=c("val"))))
dat.wd <- rbind(dat.wd,
                data.frame(deg=4, read.table("params/stat_deg4/wd_stats.txt",
                             col.names=c("val"))))
dat.wd$deg <- as.factor(dat.wd$deg)
dat.wd <- data.frame(type="Weak design", dat.wd)

dat.bitext <- data.frame(deg=2, read.table("params/stat_deg2/bitext_stats.txt",
                           col.names=c("val")))
dat.bitext <- rbind(dat.bitext,
                    data.frame(deg=3, read.table("params/stat_deg3/bitext_stats.txt",
                                 col.names=c("val"))))
dat.bitext <- rbind(dat.bitext,
                    data.frame(deg=4, read.table("params/stat_deg4/bitext_stats.txt",
                                 col.names=c("val"))))
dat.bitext$deg <- as.factor(dat.bitext$deg)
dat.bitext <- data.frame(type="Bit extractor", dat.bitext)
dat.tot <- rbind(dat.bitext, dat.wd)

# Convert the measured values from ns to us
dat.tot$val <- dat.tot$val/1000

plt <- ggplot(dat=dat.tot, aes(x=deg, y=val)) +
  geom_boxplot(aes(colour=type), outlier.size=0.5) +
  labs(y="Time[us] per component execution", x="Degree (z)", colour="Type") +
  scale_y_log10()
ggsave("paper/pictures/perf.pdf", plt)
