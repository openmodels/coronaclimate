predict.felm <- function(object, newdata, se.fit = FALSE,
                         interval = "none",
                         level = 0.95){
  if(missing(newdata)){
    stop("predict.felm requires newdata and predicts for all group effects = 0.")
  }

  tt <- terms(object)
  Terms <- delete.response(tt)
  attr(Terms, "intercept") <- 0

  m.mat <- model.matrix(Terms, data = newdata)
  m.coef <- as.numeric(object$coef)
  fit <- as.vector(m.mat %*% object$coef)
  fit <- data.frame(fit = fit)

  if(se.fit | interval != "none"){
    if(!is.null(object$clustervcv)){
      vcov_mat <- object$clustervcv
    } else if (!is.null(object$robustvcv)) {
      vcov_mat <- object$robustvcv
    } else if (!is.null(object$vcv)){
      vcov_mat <- object$vcv
    } else {
      stop("No vcv attached to felm object.")
    }
    se.fit_mat <- sqrt(diag(m.mat %*% vcov_mat %*% t(m.mat)))
  }
  if(interval == "confidence"){
    t_val <- qt((1 - level) / 2 + level, df = object$df.residual)
    fit$lwr <- fit$fit - t_val * se.fit_mat
    fit$upr <- fit$fit + t_val * se.fit_mat
  } else if (interval == "prediction"){
    stop("interval = \"prediction\" not yet implemented")
  }
  if(se.fit){
    return(list(fit=fit, se.fit=se.fit_mat))
  } else {
    return(fit)
  }
}

plot.doseresp <- function(preddf, xlabel, ylabel, indep='tasmax', mymod=NULL) {
    if (is.null(mymod))
        mymod <- mod

    plotdf <- predict.felm(mymod, preddf, interval="confidence")
    plotdf$tas <- preddf[, indep]
    ggplot(plotdf, aes(tas, fit)) +
        geom_line() + geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.5) + theme_bw() +
        scale_x_continuous(expand=c(0, 0)) + xlab(xlabel) + ylab(ylabel) +
        scale_y_continuous(labels=scales::percent)
}
