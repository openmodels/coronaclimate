outdir <- "../../results"
do.del <- F

finals <- list()
todelete <- c()

for (filename in list.files(outdir)) {
   if (any(grep("-nodel.csv", filename)))
     next
   if (any(grep("-dynamics.csv", filename)))
     next
   print(filename)
   tryCatch({
      df <- read.csv(file.path(outdir, filename))
      if (nrow(df) > 0) {
         parts <- strsplit(filename, "\\.csv")[[1]]
         if (!(parts[1] %in% names(finals)))
            finals[[parts[1]]] <- data.frame()
         finals[[parts[1]]] <- rbind(finals[[parts[1]]], df)
	 todelete <- c(todelete, file.path(outdir, filename))
      }
   }, error=function(e) {
      print("Failed!")
   })
}

for (filebase in names(finals)) {
    if (do.del)
        outpath <- file.path(outdir, paste0(filebase, ".csv"))
    else
        outpath <- file.path(outdir, paste0(filebase, "-nodel.csv"))
    write.csv(finals[[filebase]], outpath, row.names=F)
    todelete <- todelete[todelete != outpath]
}

if (do.del) {
    for (filepath in todelete) {
        file.remove(filepath)
    }
}
