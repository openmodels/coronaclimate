outdir <- "../../results"
version <- "0314"

counts <- list()
lastupds <- list()
for (filename in list.files(outdir)) {
    if (substring(filename, nchar("epimodel-") + 1, nchar("epimodel-") + nchar(version)) != version)
        next

    basename <- strsplit(filename, ".csv")[[1]][1]
    subversion <- substring(basename, nchar("epimodel-") + nchar(version) + 1, nchar(basename))
    if (!(subversion %in% names(counts))) {
        counts[[subversion]] <- 0
        lastupds[[subversion]] <- file.info(file.path(outdir, filename))$mtime
    }
    counts[[subversion]] <- counts[[subversion]] + 1
    lastupds[[subversion]] <- max(lastupds[[subversion]], file.info(file.path(outdir, filename))$mtime)
}

for (subversion in names(counts)) {
    if (counts[[subversion]] > 1)
        print(paste0(subversion, ": ", counts[[subversion]], " at ", lastupds[[subversion]]))
}
