##### Variance- Mean Dependency Plot
library(touchard)
library(ggplot2)
library(latex2exp)


r <- matrix(NA, nrow=7, ncol=11)
lambda <- seq(0.05,100,by=0.1) #values of lambda
delta <- c(-7,seq(-5,10,by=5)) #selected values of delta

datalist = list()

#calculate mu and variance for each values of lambda and delta
for ( j in delta) {
  ex = c()
  va = c()
  for ( i in lambda ) {
    t0 <- tau(i,j)
    t1 <- tau(i,j+1)
    t2 <- tau(i,j+2)
    ex <- c(ex,(t1/t0) - 1)
    va <- c(va,(t2/t0) -(t1/t0)*(t1/t0))
  }
  dat <- data.frame(ex, va)
  dat$j <- j 
  datalist[[which(j==delta)]] <- dat # add it to your list
}

big_data = do.call(rbind, datalist)
#plot the data
p <- ggplot(big_data, aes(x=ex, y=va, color=as.factor(j)))+
  geom_line(size=0.8) +
  xlab(TeX('$\\mu$')) +
  ylab(TeX('$\\sigma^2$')) +
  ggtitle(TeX(' ')) +
  coord_cartesian(ylim=c(-1, 70),xlim=c(-1, 50)) +
  guides(color=guide_legend(title=NULL)) +
  scale_color_discrete(labels=lapply(sprintf('$\\delta = %d$', delta), TeX)) +
  theme_bw()
print(p)