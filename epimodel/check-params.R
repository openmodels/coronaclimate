outdir <- "../../results"
basename <- "epimodel-0314-full2.csv"

allparams <- c()
allmus <- c()
allrhats <- c()

for (filename in list.files(outdir)) {
   if (substring(filename, 1, nchar(basename)) != basename)
      next
     
   print(filename)
   tryCatch({
      df <- read.csv(file.path(outdir, filename))
      if (nrow(df) > 0) {
         allparams <- c(allparams, as.character(df$param))
	 allmus <- c(allmus, df$mu)
	 allrhats <- c(allrhats, df$rhat)
      }
   }, error=function(e) {
      print("Failed!")
   })
}

results <- data.frame()
for (param in unique(allparams)) {
   mumu <- mean(allmus[allparams == param], na.rm=T)
   sdmu <- sd(allmus[allparams == param], na.rm=T)
   rhat <- median(allrhats[allparams == param], na.rm=T)
   results <- rbind(results, data.frame(param, mumu, sdmu, rhat))
}

print(results)
