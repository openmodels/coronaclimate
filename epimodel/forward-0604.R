## Just model the derivative
dforward <- function(ss, ee1, ee2, ii1, ii2, qq, fout, Nin, alpha, beta, sigma, gamma, kappa, N) {
    dlogbeta <- rnorm(1, 0, alpha)
    dss <- -beta*ss*(ii1 + ii2) / N
    dee1 <- (1 - fout)*beta*ss*(ii1 + ii2)/N + Nin - 2*sigma*ee1
    dee2 <- 2*sigma*ee1 - 2*sigma*ee2
    dii1 <- 2*sigma*ee2 - 2*gamma*ii1
    dii2 <- 2*gamma*ii1 - 2*gamma*ii2
    dqq <- 2*sigma*ee2*exp(-gamma*kappa) - kappa*qq
    ddd <- 2*sigma*ee2*exp(-gamma*kappa)
    dcc <- kappa*qq

    list(dlogbeta, dss, dee1, dee2, dii1, dii2, dqq, ddd, dcc)
}


alpha <- 0.395 ## beta (transmission) volatility
sigma <- 1 / 5.2 ## rate of becoming symptomatic (1 / incubation period)
gamma <- 1 / 2.9 ## rate of isolation/recovery (1/delay from onset-to-hospitalisation)
kappa <- 1 / 6.1 ## rate of reporting (1/delay from onset-to-confirmation)
N <- 11e6
ss0 <- N
ee10 <- 1
r0 <- 2.5

beta0 <- r0 / gamma ## beta(t) is transmission rate

dforward(ss0, ee10, 0, 0, 0, 0, 0, 0, alpha, beta0, sigma, gamma, kappa, N)


## Questions:
## At individual level, Erland distributed. Why at the aggregate level?
## How account for change in reporting regime?
## Explain confirmed cases expression, and why not equal to reported confirmed cases.
