setwd("~/Dropbox/Coronavirus and Climate")

##source("code/analysis/prepare.R")
load("cases/panel-prepped.RData")

growth.tree <- function(sequence) {
    varprior <- var(sequence, na.rm=T)
    priorweight <- 14 # this long for policy to take effect
    branches <- list()
    for (ss in 1:length(sequence))
        if (!is.na(sequence[ss]))
            branches[[ss]] = list(mean=sequence[ss], evar=varprior, values=sequence[ss], children=NULL, dnorm=Inf)

    while (length(branches) > 1) {
        bestmatch <- NULL
        bestdnorm <- -Inf
        for (ss in 2:length(branches)) {
            if (is.null(branches[[ss-1]]) || is.null(branches[[ss]]))
                next
            thisdnorm <- dnorm(branches[[ss-1]]$mean, branches[[ss]]$mean, branches[[ss]]$evar, log=T) +
                dnorm(branches[[ss]]$mean, branches[[ss-1]]$mean, branches[[ss-1]]$evar, log=T)
            if (thisdnorm > bestdnorm) {
                bestmatch <- ss
                bestdnorm <- thisdnorm
            }
        }

        if (is.null(bestmatch))
            break

        newbranches <- list()
        if (ss > 2) {
            for (ss in 1:(bestmatch - 2)) {
                newbranches[[ss]] <- branches[[ss]]
            }
        }
        allvalues <- c(branches[[bestmatch-1]]$values, branches[[bestmatch]]$values)
        ## Apply Scaled inverse chi-squared prior and get modal variance
        nu <- priorweight + length(allvalues)
        tausqr <- (priorweight * varprior + sum((allvalues - mean(allvalues))^2)) / (priorweight + length(allvalues))
        modevar <- nu * tausqr / (nu - 2)

        children <- list()
        children[[1]] <- branches[[bestmatch-1]]
        children[[2]] <- branches[[bestmatch]]
        newbranches[[bestmatch-1]] <- list(mean=mean(allvalues), evar=modevar, values=allvalues, children=children, dnorm=bestdnorm)
        if (bestmatch != length(branches)) {
            for (ss in (bestmatch+1):length(branches)) {
                newbranches[[ss-1]] <- branches[[ss]]
            }
        }
        branches <- newbranches
    }

    branches
}

extract.segments <- function(branches, t0, depth) {
    plotdf <- data.frame()
    for (ss in 1:length(branches)) {
        if (is.null(branches[[ss]]))
            next

        plotdf <- rbind(plotdf, data.frame(mean=branches[[ss]]$mean, evar=branches[[ss]]$evar, dnorm=branches[[ss]]$dnorm, t0=t0 + ss - 1, t1 = t0 + ss - 1 + length(branches[[ss]]$values) - 1, depth))
        if (!is.null(branches[[ss]]$children))
            plotdf <- rbind(plotdf, extract.segments(list(branches[[ss]]$children[[1]]), t0 + ss - 1, depth=depth+1),
                            extract.segments(list(branches[[ss]]$children[[2]]), t0 + ss - 1 + length(branches[[ss]]$children[[1]]$values), depth=depth+1))
    }

    plotdf
}

subdf = subset(df, regid == "United Kingdom  ")
sequence <- subdf$dlog
branches <- growth.tree(sequence)

plotdf <- extract.segments(branches, 1, 0)
plotdf$point <- plotdf$t0 == plotdf$t1
plotdf$xdepth <- plotdf$depth
plotdf$xdepth[plotdf$xdepth > 3] <- NA

ggplot(plotdf, aes(y=mean, colour=factor(xdepth))) +
    geom_segment(data=subset(plotdf, !point), aes(x=t0, xend=t1, yend=mean)) +
    geom_point(data=subset(plotdf, point), aes(x=t0))

ggplot(plotdf, aes(y=mean, colour=dnorm)) +
    geom_segment(data=subset(plotdf, !point), aes(x=t0, xend=t1, yend=mean)) +
    geom_point(data=subset(plotdf, point), aes(x=t0))

quantile(plotdf$dnorm[plotdf$dnorm < Inf])

extract.fes <- function(branches, limit) {
    fedf <- data.frame()
    for (ss in 1:length(branches)) {
        if (is.null(branches[[ss]]))
            next

        fedf <- rbind(fedf, extract.fes.helper(branches[[ss]], ss, ss, limit))
    }
    fedf
}

extract.fes.helper <- function(root, fe, name, limit) {
    if (is.null(root$children))
        return(data.frame(mean=root$mean, fe))

    name.left <- paste0(name, 'l')
    name.right <- paste0(name, 'r')

    ## Should my children have the same FE?
    if (root$dnorm < 2*limit)
        return(rbind(extract.fes.helper(root$children[[1]], name.left, name.left, limit),
                     extract.fes.helper(root$children[[2]], name.right, name.right, limit)))
    else
        return(rbind(extract.fes.helper(root$children[[1]], fe, name.left, limit),
                     extract.fes.helper(root$children[[2]], fe, name.right, limit)))
}

fes <- extract.fes(branches, dnorm(qnorm(.995), log=T))
fes$Date <- subdf$Date[!is.na(subdf$dlog)]
multis <- names(table(fes$fe)[table(fes$fe) > 1])
fes$fe2 <- fes$fe
fes$fe2[!(fes$fe %in% multis)] <- NA

ggplot(fes, aes(as.Date(Date), mean, colour=fe2)) +
    geom_point() + scale_x_date() + theme_bw() + xlab(NULL) + ylab("Delta log(Confirmed Cases)") +
    scale_colour_discrete(name="Fixed\nEffects")

df$fex <- ""
for (rid in unique(df$regid)) {
    print(rid)
    rows <- df$regid == rid
    subdf = df[rows,]
    sequence <- subdf$dlog
    if (sum(!is.na(sequence)) == 0)
        next
    branches <- growth.tree(sequence)

    fes <- extract.fes(branches, dnorm(qnorm(.995), log=T))
    multis <- names(table(fes$fe)[table(fes$fe) > 1])
    fes$fe2 <- fes$fe
    fes$fe2[!(fes$fe %in% multis)] <- NA

    subdffe <- ""
    subdffe[!is.na(subdf$dlog)] <- paste0(rid, '-', fes$fe2)
    subdffe[!is.na(subdf$dlog)][is.na(fes$fe2)] <- ""
    subdffe[is.na(subdffe)] <- ""

    df$fex[rows] <- subdffe
}

mod <- felm(dlog ~ t2m.pred + t2m2.pred + tp.pred | factor(regid) + factor(regid) : factor(week) + factor(superset) : factor(Date) | 0 | regid, data=df)
summary(mod)

mod <- felm(dlog ~ t2m.pred + t2m2.pred + tp.pred | factor(regid) + factor(fex) + factor(superset) : factor(Date) | 0 | regid, data=subset(df, fex != ''))
summary(mod)

## > length(unique(df$fex))
## [1] 4361
## > length(unique(df$fex))
## [1] 4349

write.csv(df[, c(1:10, ncol(df))], "policy/policyml.csv", row.names=F)

df2 <- df[df$fex != '',] %>% group_by(regid, fex) %>% summarize(dlog.mean=mean(dlog, na.rm=T), dlog.var=var(dlog, na.rm=T))

write.csv(df2, "policy/policyml-fexinfo.csv", row.names=F)
