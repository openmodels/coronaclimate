## setwd("~/research/coronavirus/code/seir-model")

model <- function(betas, gamma, N) {
    SS <- N-3
    II <- 3
    RR <- 0
    CC <- 3
    JJ <- 3
    for (tt in 1:length(betas)) {
        SS[tt+1] <- SS[tt] - betas[tt]*II[tt]*SS[tt] / N
        II[tt+1] <- II[tt] + betas[tt]*II[tt]*SS[tt] / N - gamma * II[tt]
        RR[tt+1] <- RR[tt] + gamma * II[tt]
        CC[tt+1] <- betas[tt]*II[tt]*SS[tt] / N
        JJ[tt+1] <- CC[tt+1] + (1 - 1 / 8.1) * JJ[tt]
    }
    list(SS=SS[-length(SS)], II=II[-length(II)], RR=RR[-length(RR)], CC=CC[-length(CC)], JJ=JJ[-length(JJ)])
}

res0.sir <- model(rep(.2, 100), .1, 1000)
res1.sir <- model(c(rep(.2, 20), .4, rep(.2, 79)), .1, 1000)

res0.sir$igrow <- c(NA, log(res0.sir$II[-1]) - log(res0.sir$II[-length(res0.sir$II)]))
res1.sir$igrow <- c(NA, log(res1.sir$II[-1]) - log(res1.sir$II[-length(res0.sir$II)]))

res0.sir$jgrow <- c(NA, log(res0.sir$JJ[-1]) - log(res0.sir$JJ[-length(res0.sir$JJ)]))
res1.sir$jgrow <- c(NA, log(res1.sir$JJ[-1]) - log(res1.sir$JJ[-length(res0.sir$JJ)]))

plot(res1.sir$igrow - res0.sir$igrow)
lines(res1.sir$jgrow - res0.sir$jgrow)

## Add E

model2 <- function(betas, alpha, gamma, N) {
    SS <- N-3
    EE <- 0
    II <- 3
    RR <- 0
    CC <- 3
    JJ <- 3
    for (tt in 1:length(betas)) {
        SS[tt+1] <- SS[tt] - betas[tt]*II[tt]*SS[tt] / N
        EE[tt+1] <- EE[tt] + betas[tt]*II[tt]*SS[tt] / N - alpha * EE[tt]
        II[tt+1] <- II[tt] + alpha * EE[tt] - gamma * II[tt]
        RR[tt+1] <- RR[tt] + gamma * II[tt]
        CC[tt+1] <- alpha * EE[tt]
        JJ[tt+1] <- CC[tt+1] + (1 - 1 / 8.1) * JJ[tt]
    }
    list(SS=SS[-length(SS)], EE=EE[-length(EE)], II=II[-length(II)], RR=RR[-length(RR)], CC=CC[-length(CC)], JJ=JJ[-length(JJ)])
}

res0 <- model2(rep(.2, 100), .2, .2, 1000)
res1 <- model2(c(rep(.2, 20), .4, rep(.2, 79)), .2, .2, 1000)

res0$igrow <- c(NA, log(res0$II[-1]) - log(res0$II[-length(res0$II)]))
res1$igrow <- c(NA, log(res1$II[-1]) - log(res1$II[-length(res0$II)]))

res0$jgrow <- c(NA, log(res0$JJ[-1]) - log(res0$JJ[-length(res0$JJ)]))
res1$jgrow <- c(NA, log(res1$JJ[-1]) - log(res1$JJ[-length(res0$JJ)]))

plot(res1$igrow - res0$igrow)
lines(res1$jgrow - res0$jgrow)

model3 <- function(betas, alpha, gamma, N) {
    SS <- N-3
    EE1 <- 0
    EE2 <- 0
    II1 <- 3
    II2 <- 0
    RR <- 0
    CC <- 3
    JJ <- 3
    for (tt in 1:length(betas)) {
        SS[tt+1] <- SS[tt] - betas[tt]*(II1[tt]+II2[tt])*SS[tt] / N
        EE1[tt+1] <- EE1[tt] + betas[tt]*(II1[tt]+II2[tt])*SS[tt] / N - (2 * alpha) * EE1[tt]
        EE2[tt+1] <- EE2[tt] + (2 * alpha) * EE1[tt] - (2 * alpha) * EE2[tt]
        II1[tt+1] <- II1[tt] + (2 * alpha) * EE2[tt] - (2 * gamma) * II1[tt]
        II2[tt+1] <- II2[tt] + (2 * gamma) * II1[tt] - (2 * gamma) * II2[tt]
        RR[tt+1] <- RR[tt] + (2 * gamma) * II2[tt]
        CC[tt+1] <- (2 * alpha) * EE2[tt]
        JJ[tt+1] <- CC[tt+1] + (1 - 1 / 8.1) * JJ[tt]
    }
    list(SS=SS[-length(SS)], EE=(EE1+EE2)[-length(EE1)], II=(II1+II2)[-length(II1)], RR=RR[-length(RR)], CC=CC[-length(CC)], JJ=JJ[-length(JJ)])
}

res0.full <- model3(rep(.2, 100), .2, .2, 1000)
res1.full <- model3(c(rep(.2, 20), .4, rep(.2, 79)), .2, .2, 1000)

res0.full$igrow <- c(NA, log(res0.full$II[-1]) - log(res0.full$II[-length(res0.full$II)]))
res1.full$igrow <- c(NA, log(res1.full$II[-1]) - log(res1.full$II[-length(res1.full$II)]))

res0.full$jgrow <- c(NA, log(res0.full$JJ[-1]) - log(res0.full$JJ[-length(res0.full$JJ)]))
res1.full$jgrow <- c(NA, log(res1.full$JJ[-1]) - log(res1.full$JJ[-length(res1.full$JJ)]))

model4 <- function(betas, alpha, gamma, theta, phi, N) {
    SS <- N-3
    EE1 <- 0
    EE2 <- 0
    II1 <- 3
    II2 <- 0
    RR <- 0
    TT <- 0
    WW <- 0
    DD <- 0
    CC <- 3
    JJ <- 3
    for (tt in 1:length(betas)) {
        SS[tt+1] <- SS[tt] - betas[tt]*(II1[tt]+II2[tt])*SS[tt] / N
        EE1[tt+1] <- EE1[tt] + betas[tt]*(II1[tt]+II2[tt])*SS[tt] / N - (2 * alpha) * EE1[tt]
        EE2[tt+1] <- EE2[tt] + (2 * alpha) * EE1[tt] - (2 * alpha) * EE2[tt]
        II1[tt+1] <- II1[tt] + (2 * alpha) * EE2[tt] - (2 * gamma) * II1[tt]
        II2[tt+1] <- II2[tt] + (2 * gamma) * II1[tt] - (2 * gamma) * II2[tt]
        RR[tt+1] <- RR[tt] + (2 * gamma) * II2[tt]

        TT[tt+1] <- TT[tt] + (2 * alpha) * EE1[tt] - theta * TT[tt]
        WW[tt+1] <- WW[tt] + theta * TT[tt] - phi * WW[tt]
        DD[tt+1] <- DD[tt] + phi * WW[tt]
        
        CC[tt+1] <- (2 * alpha) * EE2[tt]
        JJ[tt+1] <- CC[tt+1] + (1 - 1 / 8.1) * JJ[tt]
    }
    list(SS=SS[-length(SS)], EE=(EE1+EE2)[-length(EE1)], II=(II1+II2)[-length(II1)], RR=RR[-length(RR)], CC=CC[-length(CC)], JJ=JJ[-length(JJ)],
         DD=DD[-length(DD)])
}

res0.next <- model4(rep(.2, 100), .2, .2, .2, .2, 1000)
res1.next <- model4(c(rep(.2, 20), .4, rep(.2, 79)), .2, .2, .2, .2, 1000)

res0.next$igrow <- c(NA, log(res0.next$DD[-1]) - log(res0.next$DD[-length(res0.next$DD)]))
res1.next$igrow <- c(NA, log(res1.next$DD[-1]) - log(res1.next$DD[-length(res1.next$DD)]))

model2t <- function(betas, alpha, gamma, theta, phi, N) {
    SS <- N-3
    EE <- 0
    II <- 3
    RR <- 0
    TT <- 0
    WW <- 0
    DD <- 0
    CC <- 3
    JJ <- 3
    for (tt in 1:length(betas)) {
        SS[tt+1] <- SS[tt] - betas[tt]*II[tt]*SS[tt] / N
        EE[tt+1] <- EE[tt] + betas[tt]*II[tt]*SS[tt] / N - alpha * EE[tt]
        II[tt+1] <- II[tt] + alpha * EE[tt] - gamma * II[tt]
        RR[tt+1] <- RR[tt] + gamma * II[tt]

        TT[tt+1] <- TT[tt] + alpha * EE[tt] - theta * TT[tt]
        WW[tt+1] <- WW[tt] + theta * TT[tt] - phi * WW[tt]
        DD[tt+1] <- DD[tt] + phi * WW[tt]

        CC[tt+1] <- alpha * EE[tt]
        JJ[tt+1] <- CC[tt+1] + (1 - 1 / 8.1) * JJ[tt]
    }
    list(SS=SS[-length(SS)], EE=EE[-length(EE)], II=II[-length(II)], RR=RR[-length(RR)], CC=CC[-length(CC)], JJ=JJ[-length(JJ)],
         DD=DD[-length(DD)])
}

res0.t <- model2t(rep(.2, 100), .2, .2, .2, .2, 1000)
res1.t <- model2t(c(rep(.2, 20), .4, rep(.2, 79)), .2, .2, .2, .2, 1000)

res0.t$igrow <- c(NA, log(res0.t$DD[-1]) - log(res0.t$DD[-length(res0.t$DD)]))
res1.t$igrow <- c(NA, log(res1.t$DD[-1]) - log(res1.t$DD[-length(res1.t$DD)]))

model.t <- function(betas, gamma, theta, phi, N) {
    SS <- N-3
    II <- 3
    RR <- 0
    TT <- 0
    WW <- 0
    DD <- 0
    CC <- 3
    JJ <- 3
    for (tt in 1:length(betas)) {
        SS[tt+1] <- SS[tt] - betas[tt]*II[tt]*SS[tt] / N
        II[tt+1] <- II[tt] + betas[tt]*II[tt]*SS[tt] / N - gamma * II[tt]
        RR[tt+1] <- RR[tt] + gamma * II[tt]

        TT[tt+1] <- TT[tt] + betas[tt]*II[tt]*SS[tt] / N - theta * TT[tt]
        WW[tt+1] <- WW[tt] + theta * TT[tt] - phi * WW[tt]
        DD[tt+1] <- DD[tt] + phi * WW[tt]

        CC[tt+1] <- betas[tt]*II[tt]*SS[tt] / N
        JJ[tt+1] <- CC[tt+1] + (1 - 1 / 8.1) * JJ[tt]
    }
    list(SS=SS[-length(SS)], II=II[-length(II)], RR=RR[-length(RR)], CC=CC[-length(CC)], JJ=JJ[-length(JJ)], DD=DD[-length(DD)])
}

res0.sirt <- model.t(rep(.2, 100), .1, .2, .2, 1000)
res1.sirt <- model.t(c(rep(.2, 20), .4, rep(.2, 79)), .1, .2, .2, 1000)

res0.sirt$igrow <- c(NA, log(res0.sirt$DD[-1]) - log(res0.sirt$DD[-length(res0.sirt$DD)]))
res1.sirt$igrow <- c(NA, log(res1.sirt$DD[-1]) - log(res1.sirt$DD[-length(res0.sirt$DD)]))

df <- data.frame(tt=1:100, digrow.sir=res1.sir$igrow - res0.sir$igrow, digrow=res1$igrow - res0$igrow,
                 digrow.full=res1.full$igrow - res0.full$igrow, digrow.next=res1.next$igrow - res0.next$igrow,
                 digrow.t=res1.t$igrow - res0.t$igrow, digrow.sirt=res1.sirt$igrow - res0.sirt$igrow)


library(ggplot2)

ggplot(df, aes(tt)) +
    geom_hline(yintercept=0) + 
    geom_vline(xintercept=21) + geom_text(label="Beta shock", x=21.2, y=0.10, hjust="left") +
    geom_line(aes(y=digrow.sir, colour='SIR')) +
    geom_line(aes(y=digrow.sirt, colour='SIR+T')) +
    geom_line(aes(y=digrow, colour='SEIR')) +
    geom_line(aes(y=digrow.t, colour='SEIR+T')) +
    geom_line(aes(y=digrow.full, colour='SEEIIR')) +
    geom_line(aes(y=digrow.next, colour='SEEIIR+T')) +
    ylab("Growth rate in infections") +
    theme_bw() + scale_x_continuous("Days", expand=c(0, 0), limits=c(20, 60)) +
    scale_colour_discrete(name="Experiment:", breaks=c('SIR', 'SIR+T', 'SEIR', 'SEIR+T', 'SEEIIR', 'SEEIIR+T'))
ggsave("../../figures/sir-growth.pdf", width=6.5, height=4)

ggplot(df, aes(tt - 21)) +
    geom_hline(yintercept=0) + 
    geom_line(aes(y=digrow.sir, colour='SIR')) +
    geom_line(aes(y=digrow.sirt, colour='SIR+T')) +
    geom_line(aes(y=digrow.full, colour='SEIR')) +
    geom_line(aes(y=digrow.next, colour='SEIR+T')) +
    ylab("Growth rate in infections") +
    theme_bw() + scale_x_continuous("Days since shock", expand=c(0, 0), limits=c(0, 40)) +
    scale_colour_discrete(name="Experiment:", breaks=c('SIR', 'SIR+T', 'SEIR', 'SEIR+T'))
ggsave("../../figures/sir-growth.pdf", width=6.5, height=4)

library(reshape2)
df2 <- subset(melt(df, 'tt', variable.name='model'), !(model %in% c('digrow', 'digrow.t')))
df2$epimod <- "SIR"
## df2$epimod[df2$model %in% c('digrow', 'digrow.t')] <- "SEIR"
df2$epimod[df2$model %in% c('digrow.full', 'digrow.next')] <- "SEIR"
df2$testing <- F
df2$testing[df2$model %in% c('digrow.sirt', 'digrow.t', 'digrow.next')] <- T

library(dplyr)
df3 <- df2 %>% group_by(epimod) %>% summarize(maxval.not=max(value[!testing], na.rm=T),
                                              cumsum.15.not=sum(value[!testing][21:35]),
                                              cumsum.30.not=sum(value[!testing][21:50]),
                                              cumsum.60.not=sum(value[!testing][21:80]),
                                              maxval.wit=max(value[testing], na.rm=T),
                                              cumsum.15.wit=sum(value[testing][21:35]),
                                              cumsum.30.wit=sum(value[testing][21:50]),
                                              cumsum.60.wit=sum(value[testing][21:80]),
                                              cumsum.30.rat=cumsum.30.wit/cumsum.15.wit,
                                              cumsum.60.rat=cumsum.60.wit/cumsum.15.wit,
                                              maxval.rat=maxval.wit/maxval.not)
library(xtable)
print(xtable(df3, digits=3), include.rownames=F)
