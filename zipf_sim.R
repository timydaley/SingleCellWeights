y = rpois(10000, 0.1*exp(rnorm(10000, mean = -1.25^2/2, sd = 1.25)))
x = y[y > 0]
x_hist = data.frame(j = sort(unique(x), decreasing = FALSE), n_j = sapply(sort(unique(x), decreasing = FALSE), function(t) sum(x == t)))

sampleFPDD <- function(kappa, theta, n){
  counts = c(1)
  current_total_counts = 1
  while(current_total_counts < n){
    # u is uniform rv in (0, 1)
    u = runif(1)
    test_val = 0.0
    for(i in 1:length(counts)){
      test_val = test_val + (counts[i] - kappa)/(current_total_counts + theta)
      if (u <  test_val){
        counts[i] = counts[i] + 1
        break
      }
    }
    if (u > test_val){
      counts = c(counts, 1)
    }
    current_total_counts = current_total_counts + 1
  }
  return(counts)
}

x = sampleFPDD(0.96, 0.2, 1000)
x_hist = data.frame(j = sort(unique(x), decreasing = FALSE), n_j = sapply(sort(unique(x), decreasing = FALSE), function(t) sum(x == t)))

n_genes = 20000
p = 1/((1:n_genes)^0.9 + 0.25)
p = p/sum(p)
x = rpois(length(p), lambda = 500*p)
x = x[-which(x == 0)]
x_hist = data.frame(j = sort(unique(x), decreasing = FALSE), n_j = sapply(sort(unique(x), decreasing = FALSE), function(t) sum(x == t)))

x.poilogFit = poilog::poilogMLE(x)
x.ztnbfit = preseqR::preseqR.ztnb.em(x_hist)
x_hist.ztnbfit = data.frame(x = 0:max(x), fitted = sapply(0:max(x), function(u) (length(x)/(1 - dnbinom(0, mu = x.ztnbfit$mu, size = x.ztnbfit$size)))*dnbinom(u, mu = x.ztnbfit$mu, size = x.ztnbfit$size)))
x_hist.poilogFit = data.frame(x = 0:max(x), fitted = sapply(0:max(x), function(u) (length(x)/(1 - poilog::dpoilog(0, mu = x.poilogFit$par["mu"], sig = x.poilogFit$par["sig"]))*poilog::dpoilog(u, mu = x.poilogFit$par["mu"], sig = x.poilogFit$par["sig"]))))

head(x_hist.ztnbfit,5)
head(x_hist.poilogFit, 5)
head(x_hist, 4)

cols = ggthemes::fivethirtyeight_pal()(3)
ggplot(data.frame(x = x), aes(x = x)) + geom_histogram(fill = "black", binwidth = 1) + theme_bw() + geom_line(data = x_hist.poilogFit, aes(x = x, y = fitted), col = cols[1]) + geom_line(data = x_hist.ztnbfit, aes(x = x, y = fitted), col = cols[2]) + ylim(0, 400) 

ggplot(data.frame(x = x), aes(x = x)) + geom_histogram(fill = "black", binwidth = 1) + theme_bw() + geom_line(data = x_hist.poilogFit, aes(x = x, y = fitted), col = cols[1]) + geom_line(data = x_hist.ztnbfit, aes(x = x, y = fitted), col = cols[2]) + xlim(0, max(x))
