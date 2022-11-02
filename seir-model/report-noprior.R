setwd("~/research/coronavirus/code/seir-model")

ols.prior.length <- 15 # for OLS compare

source("../configs.R")

filepath <- "../../results/epimodel-meta-0314-noprior-all-nobs-nodel.csv"

pdf <- data.frame()
allmu <- data.frame()

df <- read.csv(filepath)

for (ww in weather) {
    row.e <- subset(df, Country == "" & Region == "" & param == paste0('e.', ww))

    pdf <- rbind(pdf, data.frame(channel='Transmission', weather=ww, level='Global',
                                 mu=row.e$mu, ylo=row.e$ci2.5, yhi=row.e$ci97.5,
                                 y25=row.e$ci25, y75=row.e$ci75))

    row.o <- subset(df, Country == "" & Region == "" & param == paste0('o.', ww))

    pdf <- rbind(pdf, data.frame(channel='Detection', weather=ww, level='Global',
                                 mu=row.o$mu, ylo=row.o$ci2.5, yhi=row.o$ci97.5,
                                 y25=row.o$ci25, y75=row.o$ci75))

    rows.e <- subset(df, Country != "" & Region == "" & group == "Combined" & param == paste0('e.', ww))
    rows.o <- subset(df, Country != "" & Region == "" & group == "Combined" & param == paste0('o.', ww))

    allmu <- rbind(allmu, data.frame(channel='Transmission', weather=ww, level='Country', mu=rows.e$mu))
    allmu <- rbind(allmu, data.frame(channel='Detection', weather=ww, level='Country', mu=rows.o$mu))
}

library(ggplot2)

pdf$weather <- factor(pdf$weather, levels=c('t2m', 'utci', 'ssrd', 'tp'))

ols <- data.frame(weather, mu=ols.priors.mu, sd=ols.priors.se)
ols$weather <- factor(ols$weather, levels=c('t2m', 'utci', 'ssrd', 'tp'))

gp <- ggplot(ols, aes(weather, mu)) +
    facet_grid(channel ~ ., scales='free_y') +
    geom_linerange(aes(ymin=mu - sd * 1.96, ymax=mu + sd * 1.96, colour='OLS 15 day'), position=position_nudge(x=-.3)) +
    geom_point(aes(colour='OLS 15 day'), position=position_nudge(x=-.3)) +
    geom_boxplot(data=pdf, stat="identity", aes(lower=y25, upper=y75, middle=mu, ymin=ylo, ymax=yhi,
                                                colour='Epi. global'), width=.25) +
    geom_violin(data=allmu, aes(fill='Epi. countries'), scale='width', width=.25, position=position_nudge(x=.3)) +
    theme_bw() + xlab("Weather variable") + ylab("Normalized weather response") +
    scale_colour_discrete(name=NULL, breaks=c('OLS 15 day', 'Epi. global')) +
    scale_fill_discrete(name=NULL) + theme(legend.position="bottom", legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm")) +
    guides(colour=guide_legend(order=1), fill=guide_legend(order=2))
ggsave("~/Dropbox/Coronavirus and Climate/figures/epimodel-noprior.pdf", gp, width=3.8, height=4.5)
