## ---- MarkovGen ----
MakeTrans <- function(cs, A) {
    pt <- A[cs, ]
    ns <- which(rmultinom(1, 1, pt) == 1)
    return(ns)
}


MarkovGen <- function(A, istart, N) {
    # A :transition matrix
    # istart : initial distribution
    # N : number of element of Markov chain to output

    out <- rep(0, N)
    s <- which(rmultinom(1, 1, istart) == 1)
    out[1] <- s
    for (i in 2:N) {
        s <- MakeTrans(s, A)
        out[i] <- s
    }

    return(out)
}

MarkovGen1 <- function(A, istart, N) {
    stopifnot(nrow(A) == ncol(A), ncol(A) == length(istart))
    out <- rep(0, N)
    state.space <- 1:nrow(A)
    out[1] <- sample(state.space, size=1, prob=istart)

    for (i in 2:N) {
        out[i] <- sample(state.space, size=1, prob=A[out[i-1], ])
    }

    return(out)
}

InferTransProbs <- function(out.states, s, k) {
    #Infer transition probs for state s

    ni <- out.states[which(out.states == s) + 1]
    pt <- table(ni) / length(ni)

    p <- rep(0.0, k)
    p[as.numeric(names(pt))] <- pt
    return(p)
}

InferModel <- function(out.states, k) {
    # out.states: output states from a Markov chain
    # k : number of states
    istart <- rep(0.0, k)
    istart[out.states[1]] <- 1.0

    A <- t(sapply(1:k, function(x) {InferTransProbs(out.states, x, k)}))
    colnames(A) <- 1:k
    rownames(A) <- 1:k

    return(list(A=A, istart=istart))
}

ReadModel <- function(trans.f, init.f) {
    A <- read.table(trans.f)
    istart <- read.table(init.f)

    return(list(A=as.matrix(A), istart=as.matrix(istart)))
}


MarkovGenTest <- function() {
    N <- 100
    trans.f <- "data/trans.dat"
    init.f <- "data/start.dat"
    model <- ReadModel(trans.f, init.f)
    out.states <- MarkovGen(model$A, as.vector(model$istart), N)

    return(out.states)
}

MarkovInferTest <- function() {
    out.states <- read.table("data/out_states.dat")
    model <- InferModel(out.states, 5)

    return(model)
}

plotOut <- function(out) {

    return(0)
}


