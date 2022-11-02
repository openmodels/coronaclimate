setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

num1all <- 1457
num2all <- 1811
numpair <- 652

df <- read.csv("../../results/pairwise.csv")
df2 <- data.frame(param=rep(df$param, 4), mu=c(df$mumu1, df$mumu2, df$mumu1all, df$mumu2all),
                  sd=c(df$musd1, df$musd2, df$musd1all, df$musd2all),
                  weather=rep(c(T, F, T, F), each=nrow(df)),
                  paired=rep(c(T, T, F, F), each=nrow(df)))

df3 <- data.frame(param=rep(df$param, 2),
                  mu=c((df$mumu2 - df$mumu1) / df$musd1,
                       (df$mumu2all - df$mumu1all) / df$musd1all),
                  sd=c(df$musd2 / df$musd1, df$musd2all / df$musd1all),
                  paired=rep(c(T, F), each=nrow(df)))

library(ggplot2)
ggplot(df3, aes(param, mu, colour=paired)) +
    coord_flip() +
    geom_point() + geom_errorbar(aes(ymin=mu - sd, ymax=mu + sd)) +
    geom_hline(yintercept=1) + geom_hline(yintercept=-1) +
    geom_hline(aes(colour=F, yintercept=1 / sqrt(num1all))) + geom_hline(aes(colour=F, yintercept=-1 / sqrt(num1all))) +
    geom_hline(aes(colour=T, yintercept=1 / sqrt(numpair))) + geom_hline(aes(colour=T, yintercept=-1 / sqrt(numpair)))

df4 <- data.frame(param=rep(df$param, 4),
                  relative=c((df$mumu1 - df$mumu2) / df$musd2,
                             (df$mumu1all - df$mumu2all) / df$musd2all,
                             df$musd1 / df$musd2, df$musd1all / df$musd2all),
                  paired=rep(c(T, F, T, F), each=nrow(df)),
                  variable=rep(c('mu', 'mu', 'sd', 'sd'), each=nrow(df)))

source("plotlib.R")
df4$label <- get.param.labels(df4$param)

ggplot(subset(df4, !is.na(relative) & variable == 'mu'), aes(label, relative, colour=paired)) +
    coord_flip() +
    geom_point() +
    geom_hline(yintercept=0) +
    geom_hline(aes(colour=F, yintercept=1 / sqrt(num2all))) + geom_hline(aes(colour=F, yintercept=-1 / sqrt(num2all))) +
    geom_hline(aes(colour=T, yintercept=1 / sqrt(numpair))) + geom_hline(aes(colour=T, yintercept=-1 / sqrt(numpair))) +
    xlab(NULL) + ylab("Relative to estimate without weather")

ggplot(subset(df4, !is.na(relative) & variable == 'sd'), aes(label, relative, colour=paired)) +
    coord_flip() +
    geom_point() +
    geom_hline(yintercept=1) +
    geom_hline(aes(colour=F, yintercept=1 + 1 / sqrt(num2all))) + geom_hline(aes(colour=F, yintercept=1 - 1 / sqrt(num2all))) +
    geom_hline(aes(colour=T, yintercept=1 + 1 / sqrt(numpair))) + geom_hline(aes(colour=T, yintercept=1 - 1 / sqrt(numpair))) +
    xlab(NULL) + ylab("Relative to estimate without weather")

### Bootstrapping

setwd("~/research/coronavirus/code/seir-model")

df <- read.csv("../../results/pairwise-all.csv")

results <- data.frame()
for (param in unique(df$param)) {
    print(param)
    rows <- which(df$param == param)

    mu.draws <- rep(NA, 10000)
    sd.draws <- rep(NA, 10000)
    for (ii in 1:10000) {
        rr <- sample(rows, length(rows), replace=T)
        mu1 <- rnorm(length(rows), df$mu1[rr], df$sd1[rr]) #df$mu1[rr]
        mu2 <- rnorm(length(rows), df$mu2[rr], df$sd2[rr]) #df$mu2[rr]

        mu.base <- mean(mu2, na.rm=T)
        sd.base <- mean(df$sd2[rr], na.rm=T)
        mu.rel <- mean((mu1 - mu.base) / sd.base, na.rm=T)
        sd.rel <- mean(df$sd1[rr] / sd.base, na.rm=T)

        mu.draws[ii] <- mu.rel
        sd.draws[ii] <- sd.rel
    }

    mu.quant <- quantile(mu.draws, c(.025, .975))
    sd.quant <- quantile(sd.draws, c(.025, .975))
    results <- rbind(results, data.frame(param, mu.mean=mean(mu.draws), mu.cilo=mu.quant[1], mu.cihi=mu.quant[2],
                                         sd.mean=mean(sd.draws), sd.cilo=sd.quant[1], sd.cihi=sd.quant[2]))
}
## Note: fails when gets to e.

source("plotlib.R")

results$label <- get.param.labels(results$param)

library(ggplot2)

## gp <- ggplot(results, aes(label, mu.mean)) +
##     coord_flip() + theme_bw() +
##     geom_point() + geom_errorbar(aes(ymin=mu.cilo, ymax=mu.cihi)) +
##     geom_hline(yintercept=0) +
##     xlab(NULL) + ylab("Average, relative to estimate without weather")
## ggsave("../../figures/pairwise-0314-mu.pdf", gp, width=6.5, height=2.7)

## gp <- ggplot(results, aes(label, sd.mean)) +
##     coord_flip() + theme_bw() +
##     geom_point() + geom_errorbar(aes(ymin=sd.cilo, ymax=sd.cihi)) +
##     geom_hline(yintercept=1) +
##     xlab(NULL) + ylab("Uncertainty (std. dev.), relative to estimate without weather")
## ggsave("../../figures/pairwise-0314-sd.pdf", gp, width=6.5, height=2.7)

mean(df$mu1[df$param == 'omega']) / mean(df$mu2[df$param == 'omega'])
mean(df$mu1[df$param == 'logomega']) - mean(df$mu2[df$param == 'logomega'])

mean(df$mu1[df$param == 'invkappa'])
mean(df$mu2[df$param == 'invkappa'])

mean(df$mu1[df$param == 'invkappa'] + df$mu1[df$param == 'invtheta'])
mean(df$mu2[df$param == 'invkappa'] + df$mu2[df$param == 'invtheta'])

mean(subset(df, param == "error")$mu1) # full3
mean(subset(df, param == "error")$mu2) # noweather

combodf <- data.frame(label=rep(results$label, 2), mu=c(results$mu.mean, results$sd.mean),
                      cilo=c(results$mu.cilo, results$sd.cilo), cihi=c(results$mu.cihi, results$sd.cihi),
                      var=rep(c('Average', 'Uncertainty (std. dev.)'), each=nrow(results)))

gp <- ggplot(subset(combodf, label != 'Exposed Imports'), aes(label, mu)) +
    facet_wrap(~ var, scales='free_x') +
    coord_flip() + theme_bw() +
    geom_point() + geom_errorbar(aes(ymin=cilo, ymax=cihi)) +
    geom_hline(data=data.frame(var=c('Average', 'Uncertainty (std. dev.)'), base=c(0, 1)), aes(yintercept=base)) +
    xlab(NULL) + ylab("Relative to estimate without weather")
ggsave("~/Dropbox/Coronavirus and Climate/figures/pairwise-0314.pdf", gp, width=8, height=2.7)


