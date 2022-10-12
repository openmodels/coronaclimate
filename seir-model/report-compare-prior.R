setwd("~/research/coronavirus/code/epimodel")

source("../configs.R")

filepaths <- list('noprior'="../../results/epimodel-meta-0314-noprior-all-nobs-nodel.csv",
                  'full3'="../../results/epimodel-meta-0314-full3-all-nobs-nodel.csv",
                  'full'="../../results/epimodel-meta-0314-full-all-nobs-nodel.csv")

pdf <- data.frame()
allmu <- data.frame()
for (group in names(filepaths)) {
    filepath <- filepaths[[group]]
    df <- read.csv(filepath)

    for (ww in weather) {
        row.e <- subset(df, Country == "" & Region == "" & param == paste0('e.', ww))
        
        pdf <- rbind(pdf, data.frame(channel='Transmission', weather=ww, level='Global', group,
                                     mu=row.e$mu, ylo=row.e$ci2.5, yhi=row.e$ci97.5,
                                     y25=row.e$ci25, y75=row.e$ci75))

        row.o <- subset(df, Country == "" & Region == "" & param == paste0('o.', ww))

        pdf <- rbind(pdf, data.frame(channel='Detection', weather=ww, level='Global', group,
                                     mu=row.o$mu, ylo=row.o$ci2.5, yhi=row.o$ci97.5,
                                     y25=row.o$ci25, y75=row.o$ci75))
    
        rows.e <- subset(df, Country != "" & Region == "" & group == "Combined" & param == paste0('e.', ww))
        
        pdf <- rbind(pdf, data.frame(channel='Transmission', weather=ww, level='Country', group,
                                     mu=NA, ylo=quantile(rows.e$mu, .025), yhi=quantile(rows.e$mu, .975),
                                     y25=quantile(rows.e$mu, .25), y75=quantile(rows.e$mu, .75)))
    
        rows.o <- subset(df, Country != "" & Region == "" & group == "Combined" & param == paste0('o.', ww))
        
        pdf <- rbind(pdf, data.frame(channel='Detection', weather=ww, level='Country', group,
                                     mu=NA, ylo=quantile(rows.o$mu, .025), yhi=quantile(rows.o$mu, .975),
                                     y25=quantile(rows.o$mu, .25), y75=quantile(rows.o$mu, .75)))

        allmu <- rbind(allmu, data.frame(channel='Transmission', weather=ww, level='Country', group, mu=rows.e$mu))
        allmu <- rbind(allmu, data.frame(channel='Detection', weather=ww, level='Country', group, mu=rows.o$mu))
    }
}

library(ggplot2)

pdf$group.label <- ""
pdf$group.label[pdf$group == 'noprior'] <- "No prior"
pdf$group.label[pdf$group == 'full3'] <- "30 day prior"
pdf$group.label[pdf$group == 'full'] <- "60 day prior"
pdf$group.label <- factor(pdf$group.label, levels=c("No prior", "30 day prior", "60 day prior"))

pdf$weather <- factor(pdf$weather, levels=c('t2m', 'utci', 'ssrd', 'tp'))

gp <- ggplot(subset(pdf, level == 'Country'), aes(weather, mu, colour=group.label)) +
    facet_grid(channel ~ ., scales='free_y') +
    geom_boxplot(stat="identity", aes(lower=y25, upper=y75, middle=mu, ymin=ylo, ymax=yhi), width=.75) +
    theme_bw() + xlab("Weather variable") + ylab("Normalized weather response") +
    scale_colour_discrete(name=NULL) + theme(legend.position="bottom")
    #geom_point(position=position_dodge(width=.5)) +
    #geom_linerange(aes(ymin=ylo, ymax=yhi), position=position_dodge(width=.5))
ggsave("~/Dropbox/Coronavirus and Climate/figures/epimodel-noprior.pdf", gp, width=3.8, height=4.5)

allmu$group.label <- ""
allmu$group.label[allmu$group == 'noprior'] <- "No prior"
allmu$group.label[allmu$group == 'full3'] <- "30 day prior"
allmu$group.label[allmu$group == 'full'] <- "60 day prior"
allmu$group.label <- factor(allmu$group.label, levels=c("No prior", "30 day prior", "60 day prior"))

allmu$weather <- factor(allmu$weather, levels=c('t2m', 'utci', 'ssrd', 'tp'))

gp <- ggplot(allmu, aes(weather, mu, fill=group.label)) +
    facet_grid(channel ~ ., scales='free_y') +
    geom_violin(scale='width', width=.75) +
    theme_bw() + xlab("Weather variable") + ylab("Normalized weather response") +
    scale_fill_discrete(name=NULL) + theme(legend.position="bottom")
    #geom_point(position=position_dodge(width=.5)) +
    #geom_linerange(aes(ymin=ylo, ymax=yhi), position=position_dodge(width=.5))
ggsave("~/Dropbox/Coronavirus and Climate/figures/epimodel-noprior.pdf", gp, width=3.8, height=4.5)
