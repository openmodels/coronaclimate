setwd("~/research/coronavirus/code/seir-model")

outdir <- "../../results"
textin <- "full3"

results <- data.frame()
for (filename in list.files(outdir)) {
   if (length(grep("-dynamics.csv", filename)) == 0 || length(grep(textin, filename)) == 0)
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

results$label <- "Transmission rate (beta)"
results$label[results$param == 'omega'] <- "Detection rate (omega)"

ggplot(subset(results, !invalid), aes(date, exp(logval), group=regid)) +
    facet_grid(label ~ assumps, scales="free_y") +
    geom_line(size=.1, alpha=.25) + coord_cartesian(ylim=c(1e-2, 1)) +
    scale_y_log10() + theme_bw() + scale_x_date(expand=c(0, 0))

## Combine across model runs
library(dplyr)

results.xassump <- results %>% group_by(version, regid, date, param, label) %>% summarize(allinvalid=all(invalid), logval=mean(logval[!invalid]), dowval=mean(dowval[!invalid]))

ggplot(subset(results.xassump, !allinvalid), aes(date, exp(logval), group=regid)) +
    facet_wrap(~ label, ncol=1, scales="free_y") +
    geom_line(size=.15, alpha=.25) + coord_cartesian(ylim=c(1e-2, 1)) +
    scale_y_log10() + scale_x_date(expand=c(0, 0)) + theme_bw() +
    xlab("Date in 2020") + ylab("Baseline regional rate")
ggsave("../../figures/dynamics-0314-full3.pdf", width=6.5, height=5)

## Some stats
results.stats <- results %>% group_by(version, regid, param) %>% summarize(minlogval=min(logval[!invalid], na.rm=T), maxlogval=max(logval[!invalid], na.rm=T))
results.stats$minlogval[!is.finite(results.stats$minlogval)] <- NA
results.stats$maxlogval[!is.finite(results.stats$maxlogval)] <- NA
results.stats$difflogval <- results.stats$maxlogval - results.stats$minlogval
quantile(exp(results.stats$difflogval[results.stats$param == 'omega']), na.rm=T)
mean(exp(results.stats$difflogval[results.stats$param == 'omega']) > 2, na.rm=T)
