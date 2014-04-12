DiscretiseCG <- function(f, k) {
    # k  number of bins to divide into

    cg <- as.matrix(read.table(f))
    cg <- as.vector(cg)

    ival <- (max(cg) - min(cg)) / k

    bins <- ceiling(cg / i)
    bins[bins == 0] = 1 #because ceiling(0)=0

    return(bins)
}

bins <- DiscretiseCG("data/cg.dat", 5)
write.table(bins, "data/dsCG.dat", row.names=FALSE, col.names=FALSE)

