library(ggplot2)

## Assume that 100 cases are available for detection each day
df <- data.frame(time=1:30, infected=rep(100, 30))

omega <- 0.25
df$observed0 <- df$infected * omega

## In the absence of a weather shock, growth rates are constant at 0
df$growth0 <- c(NA, log(df$observed0[2:nrow(df)]) - log(df$observed0[1:(nrow(df)-1)]))

## Now, suppose omega increases to 1 only on day 10
df$observed1 <- df$infected * c(rep(omega, 9), .5, rep(omega, 20))
df$growth1 <- c(NA, log(df$observed1[2:nrow(df)]) - log(df$observed1[1:(nrow(df)-1)]))

df$observed2 <- df$infected * c(rep(omega, 9), rep(.5, 21))
df$growth2 <- c(NA, log(df$observed2[2:nrow(df)]) - log(df$observed2[1:(nrow(df)-1)]))

ggplot(df, aes(time)) +
    geom_line(aes(y=growth0, colour='Constant detection', linetype='Constant detection')) +
    geom_line(aes(y=growth1, colour='Single period increase', linetype='Single period increase')) +
    # geom_line(aes(y=growth2, colour='Level shift increase', linetype='Level shift increase')) +
    # scale_colour_manual("Experiments", breaks=c('Constant detection', 'Single period increase', 'Level shift increase'), values=1:3) +
    # scale_linetype_manual("Experiments", breaks=c('Constant detection', 'Single period increase', 'Level shift increase'), values=c('solid', 'solid', 'dashed')) +
    scale_colour_manual("Experiments", breaks=c('Constant detection', 'Single period increase'), values=1:2) +
    scale_linetype_manual("Experiments", breaks=c('Constant detection', 'Single period increase'), values=c('solid', 'solid')) +
    theme_bw() + xlab("Time in Days") + ylab("Growth rate in observed cases") +
    scale_x_continuous(expand=c(0, 0))

## Assume that 100 new cases, but a pool for observations
df <- data.frame(time=1:30, infected=rep(100, 30))

omega <- 0.25
df$observable0 <- 250
df$observed0 <- 62.5
for (ii in 2:30) {
    df$observed0[ii] <- omega * df$observable0[ii-1]
    df$observable0[ii] <- .8 * (df$observable0[ii-1] - df$observed0[ii]) + df$infected[ii-1]
}

df$growth0 <- c(NA, log(df$observed0[2:nrow(df)]) - log(df$observed0[1:(nrow(df)-1)]))

## Now, suppose omega increases to 0.6 only on day 10
df$observable1 <- 250
df$observed1 <- 62.5
for (ii in 2:30) {
    if (ii == 10)
        df$observed1[ii] <- 0.5 * df$observable1[ii-1]
    else
        df$observed1[ii] <- omega * df$observable1[ii-1]
    df$observable1[ii] <- .8 * (df$observable1[ii-1] - df$observed1[ii]) + df$infected[ii-1]
}

df$growth1 <- c(NA, log(df$observed1[2:nrow(df)]) - log(df$observed1[1:(nrow(df)-1)]))

df$observable2 <- 250
df$observed2 <- 62.5
for (ii in 2:30) {
    if (ii >= 10)
        df$observed2[ii] <- 0.5 * df$observable2[ii-1]
    else
        df$observed2[ii] <- omega * df$observable2[ii-1]
    df$observable2[ii] <- .8 * (df$observable2[ii-1] - df$observed2[ii]) + df$infected[ii-1]
}

df$growth2 <- c(NA, log(df$observed2[2:nrow(df)]) - log(df$observed2[1:(nrow(df)-1)]))

ggplot(df, aes(time)) +
    geom_line(aes(y=growth0, colour='Constant detection', linetype='Constant detection')) +
    geom_line(aes(y=growth1, colour='Single period increase', linetype='Single period increase')) +
    # geom_line(aes(y=growth2, colour='Level shift increase', linetype='Level shift increase')) +
    # scale_colour_manual("Experiments", breaks=c('Constant detection', 'Single period increase', 'Level shift increase'), values=1:3) +
    # scale_linetype_manual("Experiments", breaks=c('Constant detection', 'Single period increase', 'Level shift increase'), values=c('solid', 'solid', 'dashed')) +
    scale_colour_manual("Experiments", breaks=c('Constant detection', 'Single period increase'), values=1:2) +
    scale_linetype_manual("Experiments", breaks=c('Constant detection', 'Single period increase'), values=c('solid', 'solid')) +
    theme_bw() + xlab("Time in Days") + ylab("Growth rate in observed cases") +
    scale_x_continuous(expand=c(0, 0))
