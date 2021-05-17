outdirs <- c("../../results")

allrhats <- list()

for (outdir in outdirs) {
for (filename in list.files(outdir)) {
   print(filename)
   tryCatch({
      df <- read.csv(file.path(outdir, filename))
      if (nrow(df) > 0) {
         parts <- strsplit(filename, "\\.csv")[[1]]
         if (!(parts[1] %in% names(allrhats)))
            allrhats[[parts[1]]] <- c()
         allrhats[[parts[1]]] <- c(allrhats[[parts[1]]], df$rhat)
      }
   }, error=function(e) {
      print("Failed!")
   })
}
}

for (filebase in names(allrhats)) {
   print(filebase)
   print(quantile(allrhats[[filebase]], c(.05, .25, .5, .75, .95), na.rm=T))
}
