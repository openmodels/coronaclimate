outdir <- "../../results"

finals <- list()

for (filename in list.files(outdir)) {
   print(filename)
   tryCatch({
      df <- read.csv(file.path(outdir, filename))
      if (nrow(df) > 0) {
         parts <- strsplit(filename, "\\.csv")[[1]]
         if (!(parts[1] %in% names(finals)))
            finals[[parts[1]]] <- data.frame()
         finals[[parts[1]]] <- rbind(finals[[parts[1]]], df)
      }
   }, error=function(e) {
      print("Failed!")
   })
}

for (filebase in names(finals)) {
   outpath <- file.path(outdir, paste0(filebase, "-nodel.csv"))
   write.csv(finals[[filebase]], outpath, row.names=F)
}
