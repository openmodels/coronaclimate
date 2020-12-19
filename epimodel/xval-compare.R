setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

library(dplyr)

df.xval <- read.csv("../../code/automated/crossval_length5_lastday_countries.csv")
df.meta <- read.csv("../../results/epimodel-meta-1111-mixed-all-nobs.csv")

subset(df.meta, Country == "Peru" & Region == "" & Locality == "" & group == "Combined" & param == "alpha")



df.meta2 <- df.meta %>% filter(Country != "" & Region == "" & Locality == "" & group == "Combined") %>% group_by(param, Country) %>% summarize(mu=mu[nobs == max(nobs)], nobs=max(nobs))
df.meta3 <- df.meta2 %>% group_by(Country) %>% summarize(delay.trans=mu[param == 'invgamma'][1] + mu[param == 'invsigma'][1], delay.detect=mu[param == 'invkappa'][1] + mu[param == 'invtheta'][1], regions=max(nobs))

df <- df.xval %>% left_join(df.meta3)

library(ggplot2)

ggplot(df, aes(LAST_MAX, delay.trans, size=regions)) +
    geom_point()

ggplot(df, aes(LAST_MAX, delay.detect)) +
    geom_point(aes(size=regions)) + geom_smooth(method='lm') + theme_bw()

ggplot(df, aes(LAST_MAX, delay.detect)) +
    geom_point(aes(size=regions)) + geom_text(aes(label=Country)) +
    geom_smooth(method='lm') + theme_bw()

cor(df$LAST_MAX, df$delay.detect)

summary(lm(LAST_MAX ~ delay.detect, data=df))

## See how our values compare

df.xval <- read.csv("../../results/crossval_results_oneperiod_countries.csv")

df <- data.frame()
for (Country in unique(df.meta3$Country)) {
    options <- subset(df.xval, country == Country & first - last + 1 == 5)
    if (all(is.na(options$rsqr)))
        next
    last <- round(df.meta3$delay.detect[df.meta3$Country == Country] - 2.5)
    score1 <- rank(options$rsqr)[options$last == last] / nrow(options)
    if (length(score1) == 0)
        next
    score2 <- options$rsqr[options$last == last] / max(options$rsqr)
    score3 <- (options$rsqr[options$last == last] - min(options$rsqr)) / (max(options$rsqr) - min(options$rsqr))
    df <- rbind(df, data.frame(Country, score1, score2, score3))
}

ggplot(df, aes(score1)) +
    geom_histogram()

mean(df$score1)

ggplot(df, aes(score2)) +
    geom_histogram()

ggplot(df, aes(score3)) +
    geom_histogram()
