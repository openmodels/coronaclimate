outdir <- "../../results"

finals <- list()

for (filename in list.files(outdir)) {
    tryCatch({
        df <- read.csv(file.path(outdir, filename))
        if (nrow(df) > 0) {
            parts <- strsplit(filename, "\\.csv")[[1]]
            if (!(parts[1] %in% names(finals)))
                finals[[parts[1]]] <- c()
            finals[[parts[1]]] <- rbind(finals[[parts[1]]], df$rsqr)
        }
    }, error=function(e) {
        print("Failed!")
    })
}

for (filebase in names(finals)) {
    print(filename)
    print(quantile(finals[[filename]]))
}
