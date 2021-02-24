forward <- function(data, params) {
### Data:

    ## int<lower=0> T; // time periods
    ## int<lower=0> N; // population
    ## int<lower=0> K; // # weather preds

    ## matrix[T, K] weather;
    ## real ii_init; // usually 0, so provide as known

### Parameters:

    ## // parameters
    ## vector[T] logbeta;
    ## vector[T-1] logomega;

    ## real<lower=2, upper=17> invsigma; // below 2 and daily step doesn't work
    ## real<lower=2, upper=17> invgamma; // below 2 and daily step doesn't work
    ## real<lower=1, upper=17> invkappa; // below 1 and daily step doesn't work
    ## real<lower=1, upper=17> invtheta; // below 1 and daily step doesn't work

    ## real<lower=0, upper=.1> deathrate; // rate of death
    ## real<lower=-.01, upper=0> deathlearning; // decreases deathrate
    ## real<lower=0, upper=1> deathomegaplus; // additional rate of reported deaths

    ## // effect of weather
    ## vector<lower=-.1, upper=.1>[K] effect;
    ## vector<lower=-.1, upper=.1>[K] omegaeffect;
    ## vector<lower=-.1, upper=.1>[6] doweffect6;
    ## vector<lower=-.1, upper=.1>[6] dowomegaeffect6;

    ## // latent variables
    ## vector<lower=0>[T-1] eein;

### Variables defined

    ## // latent variables
    ## vector<lower=0>[T] ss; // susceptible
    ## vector<lower=0>[T-1] new_ee1; // newly exposed
    ## vector<lower=0>[T] ee1; // exposed
    ## vector<lower=0>[T] ee2;
    ## vector<lower=0>[T] ii1; // infected
    ## vector<lower=0>[T] ii2;
    ## vector<lower=0>[T] qq; // waiting to be tested
    ## vector<lower=0>[T] rr; // waiting to be reported

    ## vector[T-1] omega;
    ## vector<lower=0>[T-1] dcc; // confirmed cases
    ## vector<lower=0>[T-1] ddeaths; // deaths

### Forward simulation

    doweffect <- rep(0, 7)
    dowomegaeffect <- rep(0, 7)
    for (dd in 1:6) {
        doweffect[dd] <- params$doweffect6[dd];
        dowomegaeffect[dd] <- params$dowomegaeffect6[dd];
    }

    ss <- rep(NA, data$N)
    new_ee1 <- rep(NA, data$N - 1)
    ee1 <- rep(NA, data$N)
    ee2 <- rep(NA, data$N)
    ii1 <- rep(NA, data$N)
    ii2 <- rep(NA, data$N)
    qq <- rep(NA, data$N)
    rr <- rep(NA, data$N)

    omega <- rep(NA, data$N - 1)
    dcc <- rep(NA, data$N - 1)
    ddeaths <- rep(NA, data$N - 1)

    ss[1] <- data$N;
    ee1[1] <- data$ii_init;
    ee2[1] <- data$ii_init;
    ii1[1] <- data$ii_init;
    ii2[1] <- data$ii_init;
    qq[1] <- data$ii_init;
    rr[1] <- data$ii_init;
    for (tt in 2:data$T) {
        new_ee1[tt-1] <- exp(params$logbeta[tt-1] + doweffect[1 + (tt %% 7)] + sum(data$weather[tt-1] * params$effect))*ss[tt-1]*(ii1[tt-1] + ii2[tt-1]) / data$N;
        ss[tt] <- ss[tt-1] - new_ee1[tt-1];
        ee1[tt] <- ee1[tt-1] + new_ee1[tt-1] - 2*ee1[tt-1]/params$invsigma + params$eein[tt-1];
        ee2[tt] <- ee2[tt-1] + 2*ee1[tt-1]/params$invsigma - 2*ee2[tt-1]/params$invsigma;
        ii1[tt] <- ii1[tt-1] + 2*ee2[tt-1]/params$invsigma - 2*ii1[tt-1]/params$invgamma;
        ii2[tt] <- ii2[tt-1] + 2*ii1[tt-1]/params$invgamma - 2*ii2[tt-1]/params$invgamma;

        qq[tt] <- qq[tt-1] + new_ee1[tt-1] - qq[tt-1]/params$invkappa;

        omega[tt-1] <- (exp(params$logomega[tt-1]) / (1 + exp(params$logomega[tt-1]))) * exp(dowomegaeffect[1 + (tt %% 7)] + sum(data$weather[tt-1] * params$omegaeffect));
        rr[tt] <- rr[tt-1] + omega[tt-1] * qq[tt-1]/params$invkappa - rr[tt-1]/params$invtheta;

        dcc[tt-1] <- rr[tt-1]/params$invtheta;
        ddeaths[tt-1] <- (2*ii2[tt-1]/params$invgamma) * params$deathrate * exp(tt * params$deathlearning) * (omega[tt-1] + (1 - omega[tt-1]) * params$deathomegaplus);
    }

    ## return(list(ss=ss, new_ee1=new_ee1, ee1=ee1, ee2=ee2, ii1=ii1, ii2=ii2, qq=qq, rr=rr, omega=omega, dcc=dcc, ddeaths=ddeaths))
    return(data.frame(TT=1:data$T, ss=ss, new_ee1=c(0, new_ee1), ee1=ee1, ee2=ee2, ii1=ii1, ii2=ii2, qq=qq, rr=rr, omega=c(0, omega), dcc=c(0, dcc), ddeaths=c(0, ddeaths)))
}

data <- list(T=365, N=1, K=0, weather=matrix(0, 365, 0), ii_init=1 / 10e6)
params <- list(logbeta=rep(-0.497, 365), logomega=rep(-2.767, 364), invsigma=4.624, invgamma=2.179, invkappa=3.068,
               invtheta=2.129, deathrate=0.001, deathlearning=-0.006, deathomegaplus=0.116, effect=c(), omegaeffect=c(),
               doweffect6=rep(0, 6), dowomegaeffect6=rep(0, 6), eein=rep(0, 364))

output <- forward(data, params)
output$cases <- cumsum(output$dcc)
output$deaths <- cumsum(output$ddeaths)
output$ee <- output$ee1 + output$ee2
output$ii <- output$ii1 + output$ii2

library(reshape2)

output2 <- melt(output, 'TT')

library(ggplot2)

ggplot(subset(output2, !(variable %in% c('new_ee1', 'ee1', 'ee2', 'ii1', 'ii2', 'dcc', 'ddeaths', 'ss'))), aes(TT, value, colour=variable)) +
    geom_line()

output2$label <- as.character(output2$variable)
output2$label[output2$label == 'ee'] <- "Exposed"
output2$label[output2$label == 'ii'] <- "Infected"
output2$label[output2$label == 'qq'] <- "Untested"
output2$label[output2$label == 'rr'] <- "Unreported"
output2$label[output2$label == 'cases'] <- "Total Reported"
output2$label <- factor(output2$label, levels=c("Exposed", "Infected", "Untested", "Unreported", "Total Reported"))

ggplot(subset(output2, !(variable %in% c('new_ee1', 'ee1', 'ee2', 'ii1', 'ii2', 'dcc', 'ddeaths', 'ss', 'deaths', 'omega'))), aes(TT, value, colour=label)) +
    geom_line() + scale_y_continuous(labels=scales::percent) + xlab("Day since start of disease (1 infected per 10 million)") +
    theme_bw() + scale_colour_discrete(name=NULL) + ylab("Percent of population")

weather <- matrix(0, 365, 5)
weather[1, ] <- rnorm(5)
for (ii in 2:365)
    weather[ii, ] <- .9 * weather[ii - 1, ] + .1 * rnorm(5)
effect <- c(0.028, -0.016, 0.007, 0.008, 0.007)
omegaeffect <- c(0.009, -0.030, -0.010, -0.012, -0.013)

data <- list(T=365, N=1, K=0, weather=weather, ii_init=1 / 10e6)
params <- list(logbeta=rep(-0.497, 365), logomega=rep(-2.767, 364), invsigma=4.624, invgamma=2.179, invkappa=3.068,
               invtheta=2.129, deathrate=0.001, deathlearning=-0.006, deathomegaplus=0.116, effect=effect,
               omegaeffect=omegaeffect, doweffect6=rep(0, 6), dowomegaeffect6=rep(0, 6), eein=rep(0, 364))

output <- forward(data, params)
output$cases <- cumsum(output$dcc)
output$deaths <- cumsum(output$ddeaths)
output$ee <- output$ee1 + output$ee2
output$ii <- output$ii1 + output$ii2

output2 <- melt(output, 'TT')

output2$label <- as.character(output2$variable)
output2$label[output2$label == 'ee'] <- "Exposed"
output2$label[output2$label == 'ii'] <- "Infected"
output2$label[output2$label == 'qq'] <- "Untested"
output2$label[output2$label == 'rr'] <- "Unreported"
output2$label[output2$label == 'cases'] <- "Total Reported"
output2$label <- factor(output2$label, levels=c("Exposed", "Infected", "Untested", "Unreported", "Total Reported"))

ggplot(subset(output2, !(variable %in% c('new_ee1', 'ee1', 'ee2', 'ii1', 'ii2', 'dcc', 'ddeaths', 'ss', 'deaths', 'omega'))), aes(TT, value, colour=label)) +
    geom_line() + scale_y_continuous(labels=scales::percent) + xlab("Day since start of disease (1 infected per 10 million)") +
    theme_bw() + scale_colour_discrete(name=NULL) + ylab("Percent of population")

## Determine effect of pulses

get.dlog <- function(data, params) {
    outputs <- forward(data, params)
    synthetic.base <- rep(0, 365)
    for (ii in 2:365)
        synthetic.base[ii] <- outputs$dcc[ii-1] + (1 - 1 / 8.1) * synthetic.base[ii-1]
    diff(log(synthetic.base))
}

get.impresp <- function(baseline, pulsetime, params) {
    weather <- matrix(0, 365, 1)
    weather[pulsetime, 1] <- 1

    data <- list(T=365, N=1, K=1, weather=weather, ii_init=1 / 10e6)

    perturbed <- get.dlog(data, params)
    perturbed - baseline
}

weather <- matrix(0, 365, 1)
data <- list(T=365, N=1, K=1, weather=weather, ii_init=1 / 10e6)
params <- list(logbeta=rep(-0.497, 365), logomega=rep(-2.767, 364), invsigma=4.624, invgamma=2.179, invkappa=3.068,
               invtheta=2.129, deathrate=0.001, deathlearning=-0.006, deathomegaplus=0.116, effect=1,
               omegaeffect=1, doweffect6=rep(0, 6), dowomegaeffect6=rep(0, 6), eein=rep(0, 364))

baseline <- get.dlog(data, params)

results <- data.frame()
results.timeinv <- data.frame()
for (effect in c('beta', 'omega', 'both')) {
    for (pulsetime in c(100, 200, 300)) {
        if (effect == 'beta')
            params <- list(logbeta=rep(-0.497, 365), logomega=rep(-2.767, 364), invsigma=4.624, invgamma=2.179, invkappa=3.068,
                           invtheta=2.129, deathrate=0.001, deathlearning=-0.006, deathomegaplus=0.116, effect=1,
                           omegaeffect=0, doweffect6=rep(0, 6), dowomegaeffect6=rep(0, 6), eein=rep(0, 364))
        else if (effect == 'omega')
            params <- list(logbeta=rep(-0.497, 365), logomega=rep(-2.767, 364), invsigma=4.624, invgamma=2.179, invkappa=3.068,
                           invtheta=2.129, deathrate=0.001, deathlearning=-0.006, deathomegaplus=0.116, effect=0,
                           omegaeffect=1, doweffect6=rep(0, 6), dowomegaeffect6=rep(0, 6), eein=rep(0, 364))
        else
            params <- list(logbeta=rep(-0.497, 365), logomega=rep(-2.767, 364), invsigma=4.624, invgamma=2.179, invkappa=3.068,
                           invtheta=2.129, deathrate=0.001, deathlearning=-0.006, deathomegaplus=0.116, effect=1,
                           omegaeffect=1, doweffect6=rep(0, 6), dowomegaeffect6=rep(0, 6), eein=rep(0, 364))

        impresp <- get.impresp(baseline, pulsetime, params)

        results <- rbind(results, data.frame(effect, pulsetime, time=1:364, impresp=impresp))
        results.timeinv <- rbind(results.timeinv, data.frame(effect, pulsetime, time=1:30, impresp=impresp[pulsetime:(pulsetime+29)]))
    }
}

results$effect <- factor(results$effect, levels=c('beta', 'omega', 'both'))

ggplot(results, aes(time, impresp)) +
    facet_grid(effect ~ pulsetime) +
    geom_line() + theme_bw()

results.timeinv$effect <- factor(results.timeinv$effect, levels=c('beta', 'omega', 'both'))

ggplot(results.timeinv, aes(time, impresp, colour=factor(pulsetime))) +
    facet_grid(effect ~ .) +
    geom_line() + theme_bw()

## Determine imp. resp. for many different delays

results <- data.frame()
for (pulsetime in c(100, 250)) {
    for (delaysum in seq(5, 30)) {
        params <- list(logbeta=rep(-0.497, 365), logomega=rep(-2.767, 364), invsigma=delaysum * 5.2 / 8.1, invgamma=delaysum * 2.9 / 8.1, invkappa=3.068,
                       invtheta=2.129, deathrate=0.001, deathlearning=-0.006, deathomegaplus=0.116, effect=1,
                       omegaeffect=0, doweffect6=rep(0, 6), dowomegaeffect6=rep(0, 6), eein=rep(0, 364))

        weather <- matrix(0, 365, 1)
        data <- list(T=365, N=1, K=1, weather=weather, ii_init=1 / 10e6)

        baseline <- get.dlog(data, params)

        impresp <- get.impresp(baseline, pulsetime, params)
        results <- rbind(results, data.frame(delaysum, pulsetime, effect='beta', time=1:364, impresp=impresp))

        params$invsigma <- 5.2
        params$invgamma <- 2.9
        params$invkappa <- delaysum * 3.068 / (3.068 + 2.129)
        params$invtheta <- delaysum * 2.129 / (3.068 + 2.129)
        params$effect <- 0
        params$omegaeffect <- 1

        weather <- matrix(0, 365, 1)
        data <- list(T=365, N=1, K=1, weather=weather, ii_init=1 / 10e6)

        baseline <- get.dlog(data, params)

        impresp <- get.impresp(baseline, pulsetime, params)
        results <- rbind(results, data.frame(delaysum, pulsetime, effect='omega', time=1:364, impresp=impresp))
    }
}

library(ggplot2)

ggplot(subset(results, time >= 100 & pulsetime == 100), aes(time - 100, impresp, colour=delaysum, group=delaysum)) +
    facet_grid(effect ~ ., scale="free_y") +
    geom_line() + scale_x_continuous(name="Time from impulse", expand=c(0, 0)) + theme_bw()

ggplot(subset(results, time >= 250 & pulsetime == 250 & delaysum > 7), aes(time - 250, impresp, colour=delaysum, group=delaysum)) +
    facet_grid(effect ~ ., scale="free_y") +
    geom_line() + scale_x_continuous(name="Time from impulse", expand=c(0, 0)) + theme_bw()

ggplot(subset(results, time >= 100 & time <= 110 & pulsetime == 100), aes(time - 100, impresp, colour=delaysum, group=delaysum)) +
    facet_grid(effect ~ ., scale="free_y") + geom_rect(xmin=2, xmax=5, ymin=-.1, ymax=.2, alpha=.5) +
    geom_line() + scale_x_continuous(name="Time from impulse", expand=c(0, 0)) + theme_bw()

