setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

mobrecorded <- read.csv("../../results/epimodel-meta-1217-mobile-pop.csv")

compare <- read.csv("literature.csv")
compare <- subset(compare, Study != "Us" & Population != "Drop")

subset(mobrecorded, group == "Combined" & Country == "" & Region == "")

pooled <- data.frame(Study="Present study", Parameter=c("Incubation period (days)", "Reproduction number", "Fatality rate per 1000"),
                     Estimate=c(subset(mobrecorded, param == "invsigma" & group == "Combined" & Country == "" & Region == "")$mu,
                                subset(mobrecorded, param == "invgamma" & group == "Combined" & Country == "" & Region == "")$mu * exp(subset(mobrecorded, param == "logbeta" & group == "Combined" & Country == "" & Region == "")$mu),
                                subset(mobrecorded, param == "deathrate" & group == "Combined" & Country == "" & Region == "")$mu * 1000),
                     X95.CI.Low=c(subset(mobrecorded, param == "invsigma" & group == "Combined" & Country == "" & Region == "")$ci2.5,
                                  subset(mobrecorded, param == "invgamma" & group == "Combined" & Country == "" & Region == "")$ci2.5 * exp(subset(mobrecorded, param == "logbeta" & group == "Combined" & Country == "" & Region == "")$ci2.5),
                                  subset(mobrecorded, param == "deathrate" & group == "Combined" & Country == "" & Region == "")$ci2.5 * 1000),
                     X95.CI.High=c(subset(mobrecorded, param == "invsigma" & group == "Combined" & Country == "" & Region == "")$ci97.5,
                                   subset(mobrecorded, param == "invgamma" & group == "Combined" & Country == "" & Region == "")$ci97.5 * exp(subset(mobrecorded, param == "logbeta" & group == "Combined" & Country == "" & Region == "")$ci97.5),
                                   subset(mobrecorded, param == "deathrate" & group == "Combined" & Country == "" & Region == "")$ci97.5 * 1000),
                     Population='Pooled')

region.invsigma <- data.frame(Study="Present study", Parameter="Incubation period (days)",
                              Estimate=subset(mobrecorded, param == "invsigma" & group == "Combined" & Country != "" & nobs > 1 & Region == "")$mu,
                              X95.CI.Low=subset(mobrecorded, param == "invsigma" & group == "Combined" & Country != "" & nobs > 1 & Region == "")$ci2.5,
                              X95.CI.High=subset(mobrecorded, param == "invsigma" & group == "Combined" & Country != "" & nobs > 1 & Region == "")$ci97.5, Population='Region')
region.r0 <- data.frame(Study="Present study", Parameter="Reproduction number",
                        Estimate=subset(mobrecorded, param == "invgamma" & group == "Combined" & Country != "" & nobs > 1 & Region == "")$mu * exp(subset(mobrecorded, param == "logbeta" & group == "Combined" & Country != "" & nobs > 1 & Region == "")$mu),
                        X95.CI.Low=subset(mobrecorded, param == "invgamma" & group == "Combined" & Country != "" & nobs > 1 & Region == "")$ci2.5 * exp(subset(mobrecorded, param == "logbeta" & group == "Combined" & Country != "" & nobs > 1 & Region == "")$ci2.5),
                        X95.CI.High=subset(mobrecorded, param == "invgamma" & group == "Combined" & Country != "" & nobs > 1 & Region == "")$ci97.5 * exp(subset(mobrecorded, param == "logbeta" & group == "Combined" & Country != "" & nobs > 1 & Region == "")$ci97.5), Population='Region')
region.fatality <- data.frame(Study="Present study", Parameter="Fatality rate per 1000",
                              Estimate=subset(mobrecorded, param == "deathrate" & group == "Combined" & Country != "" & nobs > 1 & Region == "")$mu * 1000,
                              X95.CI.Low=subset(mobrecorded, param == "deathrate" & group == "Combined" & Country != "" & nobs > 1 & Region == "")$ci2.5 * 1000,
                              X95.CI.High=subset(mobrecorded, param == "deathrate" & group == "Combined" & Country != "" & nobs > 1 & Region == "")$ci97.5 * 1000, Population='Region')

plotdf <- rbind(compare, pooled, region.invsigma, region.r0, region.fatality)

library(ggplot2)

ggplot(plotdf, aes(Parameter, Estimate, colour=Study, shape=Population)) +
    coord_flip() +
    geom_point() +
    geom_errorbar(data=subset(plotdf, Population == 'Pooled'), aes(ymin=X95.CI.Low, ymax=X95.CI.High), position=position_dodge2(padding=.5), width=1) +
    theme_bw() + xlab(NULL) + theme(legend.box="horizontal")
