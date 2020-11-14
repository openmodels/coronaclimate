## setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

filepath <- '../../results-mixed/epimodel-1111-mixed.csv'
df <- read.csv(filepath)

combine <- function(subdf) {
    if (nrow(subdf) <= 1)
        return(subdf)

    mu <- mean(subdf$mu)
    new.sd <- sd(rnorm(1e6, subdf$mu, subdf$sd))
    ci2.5 <- median(subdf$ci2.5)
    ci25 <- median(subdf$ci25)
    ci50 <- median(subdf$ci50)
    ci75 <- median(subdf$ci75)
    ci97.5 <- median(subdf$ci97.5)
    rhat <- new.sd / sqrt(mean((subdf$sd / subdf$rhat)^2))

    return(data.frame(regid=subdf$regid[1], param=subdf$param[1], mu, sd=new.sd, ci2.5, ci25, ci50, ci75, ci97.5, rhat))
}

results <- data.frame()
for (param in unique(df$param)) {
    print(param)
    for (regid in unique(df$regid)) {
        subdf <- df[df$regid == regid & df$param == param,]
        results <- rbind(results, combine(subdf))
    }
}

write.csv(results, gsub('\\.csv', '-combo.csv', filepath), row.names=F)
