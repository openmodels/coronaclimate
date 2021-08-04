## head
RESULTPATH <- '~/Dropbox/Coronavirus and Climate/results/'
FIGUREPATH <- '~/Dropbox/Coronavirus and Climate/figures/'
CASESPATH <- '~/Dropbox/Coronavirus and Climate/cases/'
MOBILITYPATH <- '~/Dropbox/Coronavirus and Climate/mobility/'
## body
weather.delayA <- c(4, 8)
weather.delayB <- c(11, 19)
weather.subset <- c("q", "q2", "r", "r2", "ssrd", "ssrd2", "t2m", "t2m2", "tp", "tp2", "utci", "utci2")
weather.subset.stagetwo <- c("q", "q2", "r", "r2", "ssrd", "ssrd2", "t2m", "t2m2", "tp", "tp2", "utci", "utci2")
weather.delay.file <- 'crossval_length5_lastday.csv'
weather.delay.length <- 5
GAMMA <- 0.1234567901234568

weather <- c('t2m', 'tp', 'ssrd', 'utci')

if (exists("ols.prior.length")) {
    if (ols.prior.length == 15) {
        ols.priors.mu <- c(-0.0037711470543045, -0.0115532024749602, 0.0038254277658579, -0.00576977001613532)
        ols.priors.se <- c(0.00411814547761493, 0.0033460420699574, 0.00448413550578797, 0.00438879716928544)
    } else if (ols.prior.length == 30) {
        ols.priors.mu <- c(0.02773106208, -0.03134987114, -0.01027632136, -0.02558036184)
        ols.priors.se <- c(0.01222135339, 0.008463339282, 0.01140204532, 0.009860867169)
    } else if (ols.prior.length == 60) {
        ols.priors.mu <- c(0.05075947952, 0.03104620013, -0.06615251605, -0.08119237222)
        ols.priors.se <- c(0.04050863843, 0.02348355048, 0.03277420295, 0.02246738476)
    }
}
