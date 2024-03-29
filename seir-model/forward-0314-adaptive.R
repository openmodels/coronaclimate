forward.adaptive <- function(data, params, dcc.true) {
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

    ss <- rep(NA, data$T)
    new_ee1 <- rep(NA, data$T - 1)
    ee1 <- rep(NA, data$T)
    ee2 <- rep(NA, data$T)
    ii1 <- rep(NA, data$T)
    ii2 <- rep(NA, data$T)
    qq <- rep(NA, data$T)
    rr <- rep(NA, data$T)

    omega <- rep(NA, data$T - 1)
    dcc <- rep(NA, data$T - 1)
    ddeaths <- rep(NA, data$T - 1)

    extraeein <- rep(NA, data$T - 1)
    
    ss[1] <- data$N;
    ee1[1] <- data$ii_init;
    ee2[1] <- data$ii_init;
    ii1[1] <- data$ii_init;
    ii2[1] <- data$ii_init;
    qq[1] <- data$ii_init;
    rr[1] <- data$ii_init;
    for (tt in 2:data$T) {
        new_ee1[tt-1] <- exp(params$logbeta[tt-1] + doweffect[1 + (tt %% 7)] + sum(data$weather[tt-1,] * params$effect))*ss[tt-1]*(ii1[tt-1] + ii2[tt-1]) / data$N;
        ss[tt] <- ss[tt-1] - new_ee1[tt-1];
        ee1[tt] <- ee1[tt-1] + new_ee1[tt-1] - 2*ee1[tt-1]/params$invsigma + params$eein[tt-1];
        ee2[tt] <- ee2[tt-1] + 2*ee1[tt-1]/params$invsigma - 2*ee2[tt-1]/params$invsigma;
        ii1[tt] <- ii1[tt-1] + 2*ee2[tt-1]/params$invsigma - 2*ii1[tt-1]/params$invgamma;
        ii2[tt] <- ii2[tt-1] + 2*ii1[tt-1]/params$invgamma - 2*ii2[tt-1]/params$invgamma;

        qq[tt] <- qq[tt-1] + 2*ee1[tt-1]/params$invsigma - qq[tt-1]/params$invkappa;

        omega[tt-1] <- exp(params$logomega[tt-1] + dowomegaeffect[1 + (tt %% 7)] + sum(data$weather[tt-1,] * params$omegaeffect))
        rr[tt] <- rr[tt-1] + omega[tt-1] * qq[tt-1]/params$invkappa - rr[tt-1]/params$invtheta;

        dcc[tt-1] <- rr[tt-1]/params$invtheta;
        ddeaths[tt-1] <- (2*ii2[tt-1]/params$invgamma) * params$deathrate * exp(tt * params$deathlearning) * (omega[tt-1] + (1 - omega[tt-1]) * params$deathomegaplus);

        if (is.na(dcc.true[tt-1]) || round(dcc[tt-1]) == round(dcc.true[tt-1])) {
            extraeein[tt-1] <- 0
        } else if (dcc[tt-1] < dcc.true[tt-1]) {
            extraeein[tt-1] <- (dcc.true[tt-1] - dcc[tt-1]) / omega[tt-1]
        } else if (dcc[tt-1] > dcc.true[tt-1]) {
            extraeein[tt-1] <- max((dcc.true[tt-1] - dcc[tt-1]) / omega[tt-1], min(0, -ee1[tt]))
        }

        ## Correct ee1
        ee1[tt] <- ee1[tt-1] + new_ee1[tt-1] - 2*ee1[tt-1]/params$invsigma + params$eein[tt-1] + extraeein[tt-1]
    }

    return(data.frame(TT=1:data$T, ss=ss, new_ee1=c(0, new_ee1), ee1=ee1, ee2=ee2, ii1=ii1, ii2=ii2, qq=qq, rr=rr, omega=c(0, omega), dcc=c(0, dcc), ddeaths=c(0, ddeaths), extraeein=c(extraeein, NA)))
}
