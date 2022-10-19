setwd("~/research/coronavirus/code/epimodel")

outdir <- "../../results-saved"

results <- data.frame()
for (filename in list.files(outdir)) {
   if (length(grep("-dynamics.csv", filename)) == 0)
     next

   parts <- strsplit(filename, "-")[[1]]
   version <- parts[2]
   assumps <- parts[3]
   regid <- parts[5]

   df <- read.csv(file.path(outdir, filename))
   invalid <- any(df$logbeta > 0)
   date <- as.Date("2020-01-01") + (1:nrow(df)) - 1
   
   results <- rbind(results, data.frame(version, assumps, regid, invalid, date, param=rep(c('beta', 'omega'), each=nrow(df)), logval=c(df$logbeta, df$logomega), dowval=c(df$doweffect, df$dowomegaeffect)))
}

library(ggplot2)

ggplot(subset(results, !invalid), aes(date, exp(logval), group=regid)) +
    facet_grid(param ~ assumps, scales="free_y") +
    geom_line(size=.1, alpha=.25) + coord_cartesian(ylim=c(1e-3, 1)) +
    scale_y_log10()

## Combine across model runs
library(dplyr)

results.xassump <- results %>% group_by(version, regid, date, param) %>% summarize(allinvalid=all(invalid), logval=mean(logval[!invalid]), dowval=mean(dowval[!invalid]))

ggplot(subset(results.xassump, !allinvalid), aes(date, exp(logval), group=regid)) +
    facet_wrap(~ param, ncol=1, scales="free_y") +
    geom_line(size=.15, alpha=.25) + coord_cartesian(ylim=c(1e-2, 1)) +
    scale_y_log10() + scale_x_date(expand=c(0, 0)) + theme_bw() +
    xlab(NULL) + ylab("Rate")
