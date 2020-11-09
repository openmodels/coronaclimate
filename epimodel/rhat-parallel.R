outdir <- "../../results"

finals <- list()

for (filename in list.files(outdir)) {
    tryCatch({
        df <- read.csv(file.path(outdir, filename))
        if (nrow(df) > 0) {
            parts <- strsplit(filename, "\\.csv")[[1]]
            if (!(parts[1] %in% names(finals)))
                finals[[parts[1]]] <- c()
            finals[[parts[1]]] <- c(finals[[parts[1]]], df$rhat)
        }
    }, error=function(e) {
        print("Failed!")
    })
}

for (filebase in names(finals)) {
    print(filebase)
    print(quantile(finals[[filebase]], na.rm=T))
}
