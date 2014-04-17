## ---- MarkovGen ----
MarkovChain <- function(A, istart, N) {
    # Generates N elements of the Markov chain
    # specified by [A, istart]
    #
    # A transition matrix
    # istart initial distribution
    # N length of chain

    stopifnot(nrow(A) == ncol(A), ncol(A) == length(istart))
    out <- rep(0, N)
    state.space <- 1:nrow(A)
    out[1] <- sample(state.space, size=1, prob=istart)

    for (i in 2:N) {
        out[i] <- sample(state.space, size=1, prob=A[out[i-1], ])
    }

    return(out)
}

ReadModel <- function(trans.f, init.f) {
    A <- read.table(trans.f)
    istart <- read.table(init.f)

    return(list(A=as.matrix(A), istart=as.matrix(istart)))
}

MarkovGen <- function(trans.f, init.f, N) {
    # trans.f  file containing transition matrix A
    # init.f   file containing initial distribution istart

    model <- ReadModel(trans.f, init.f)
    out <- MarkovChain(model$A, model$istart, N)

    return(out)
}

# ---- Inference ----
InferTransProbs <- function(out.states, s, k) {
    # Infers transition probs for state s

    ni <- out.states[which(out.states == s) + 1]
    pt <- table(ni) / length(ni)

    p <- rep(0.0, k)
    p[as.numeric(names(pt))] <- pt
    return(p)
}

InferModel <- function(states.f, k) {
    # states.f   file containing output states from a Markov chain
    # k          number of states
    #
    # Returns:
    # [A, istart]  Inferred Markov Chain model

    out.states <- read.table("data/out_states.dat")
    istart <- rep(0.0, k)
    istart[out.states[1]] <- 1.0

    A <- t(sapply(1:k, function(s) {InferTransProbs(out.states, s, k)}))
    colnames(A) <- 1:k
    rownames(A) <- 1:k

    return(list(A=A, istart=istart))
}

#----- Tests -----
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


