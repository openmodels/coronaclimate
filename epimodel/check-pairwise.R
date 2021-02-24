outdir <- "../../results"
onebase <- "epimodel-0105noprior.csv"
twobase <- "epimodel-0105now.csv"

allparams1 <- c()
allmus1 <- c()
allsds1 <- c()
allparams1all <- c()
allmus1all <- c()
allsds1all <- c()

allparams2 <- c()
allmus2 <- c()
allsds2 <- c()
allparams2all <- c()
allmus2all <- c()
allsds2all <- c()

numpair <- 0
num1all <- 0
num2all <- 0

allfiles <- list.files(outdir)

for (filename in allfiles) {
   if (substring(filename, 1, nchar(onebase)) != onebase)
      next

   tryCatch({
      df1 <- read.csv(file.path(outdir, filename))
      if (nrow(df1) > 0) {
         allparams1all <- c(allparams1all, as.character(df1$param))
	 allmus1all <- c(allmus1all, df1$mu)
	 allsds1all <- c(allsds1all, df1$sd)
	 num1all <- num1all + 1
      }
   })
}

for (filename in allfiles) {
   if (substring(filename, 1, nchar(twobase)) != twobase)
      next

   tryCatch({
      df2 <- read.csv(file.path(outdir, filename))
      if (nrow(df2) > 0) {
         allparams2all <- c(allparams2all, as.character(df2$param))
	 allmus2all <- c(allmus2all, df2$mu)
	 allsds2all <- c(allsds2all, df2$sd)
	 num2all <- num2all + 1
      }
   })
}

for (filename in allfiles) {
   if (substring(filename, 1, nchar(onebase)) != onebase)
      next

   compname <- paste0(twobase, substring(filename, nchar(onebase) + 1, nchar(filename)))
   if (!(compname %in% allfiles))
      next
     
   print(filename)
   tryCatch({
      df1 <- read.csv(file.path(outdir, filename))
      df2 <- read.csv(file.path(outdir, compname))
      if (nrow(df1) > 0 && nrow(df2) > 0) {
         allparams1 <- c(allparams1, as.character(df1$param))
	 allmus1 <- c(allmus1, df1$mu)
	 allsds1 <- c(allsds1, df1$sd)

         allparams2 <- c(allparams2, as.character(df2$param))
	 allmus2 <- c(allmus2, df2$mu)
	 allsds2 <- c(allsds2, df2$sd)
         numpair <- numpair + 1
      }
}, error=function(e) {
      print("Failed!")
   })
}

print(c(num1all, num2all, numpair))

results <- data.frame()
for (param in unique(allparams1)) {
   mumu1 <- mean(allmus1[allparams1 == param], na.rm=T)
   musd1 <- mean(allsds1[allparams1 == param], na.rm=T)
   mumu2 <- mean(allmus2[allparams2 == param], na.rm=T)
   musd2 <- mean(allsds2[allparams2 == param], na.rm=T)

   mumu1all <- mean(allmus1all[allparams1all == param], na.rm=T)
   musd1all <- mean(allsds1all[allparams1all == param], na.rm=T)
   mumu2all <- mean(allmus2all[allparams2all == param], na.rm=T)
   musd2all <- mean(allsds2all[allparams2all == param], na.rm=T)

   results <- rbind(results, data.frame(param, mumu1, musd1, mumu2, musd2, mumu1all, musd1all, mumu2all, musd2all))
}

print(results)

write.csv(results, file.path(outdir, "pairwise.csv"), row.names=F)

write.csv(data.frame(param=allparams2, mu1=allmus1[allparams1 %in% unique(allparams2)], sd1=allsds1[allparams1 %in% unique(allparams2)], mu2=allmus2, sd2=allsds2), file.path(outdir, "pairwise-all.csv"), row.names=F)
