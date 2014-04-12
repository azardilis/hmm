MakeTrans <- function(cs, A) {
    pt <- A[cs, ]
    ns <- which(rmultinom(1, 1, pt) == 1)
    return(ns)
}

EmitSymbol <- function(cs, B) {
    pt <- B[cs, ]
    y <- which(rmultinom(1, 1, pt) == 1)

    return(y)
}

HMMGen <- function(A, B, istart, N) {
    # A  transition matrix
    # B  emission probabilities
    # istart  initial distribution
    # N  number of element of Markov chain to output

    out.es <- rep(0, N)
    out.hs <- rep(0, N)
    s <- which(rmultinom(1, 1, istart) == 1)
    out.hs[1] <- s
    out.es[1] <- which(rmultinom(1, 1, B[s, ]) == 1)

    for(i in 2:N) {
        s <- MakeTrans(s, A)
        out.hs[i] <- s
        out.es[i] <- EmitSymbol(s, B)
    }

    return(list(out.hs=out.hs, out.es=out.es))
}

ReadModel <- function(trans.f, init.f, emit.f) {
    A <- read.table(trans.f)
    istart <- read.table(init.f)
    B <- read.table(emit.f)

    return(list(A=as.matrix(A), istart=as.matrix(istart),
                B=as.matrix(B)))
}

TestHMM <- function() {
    N <- 115
    trans.f <- "data/trans_hmm.dat"
    init.f <- "data/start_hmm.dat"
    emit.f <- "data/emit.dat"
    model <- ReadModel(trans.f, init.f, emit.f)

    out <- HMMGen(model$A, model$B, model$istart, N)
    return(out)
}

PlotOut <- function(out) {
    # expect a list out containing sequences
    # of hidden and emitted states from an HMM
    # run
    T <- 1:length(out$out.hs)
    states <- as.matrix(cbind(T, out$out.hs, out$out.es))

    ltyp <- rep("solid", ncol(states)-1)
    matplot(states[, 1], states[, 2:ncol(states)], type='l', ylab='State index',
            xlab = 'time', lwd=2, lty = ltyp, ylim = c(1, 6))

    legend("topright", c("Hidden States", "Emitted States"),
           col = 1:(ncol(states)-1), lty=ltyp)
}



UpdateForward <- function(f, A, B, s, yi) {
    # Updates the forward probability for a state s
    #
    # f  forward probs at previous step
    # A  transition matrix
    # B  emission probs
    # s  state that forward is calculated for
    # y  current emitted symbol

    nf <- B[s, yi] * sum(f * A[, s])

    return(nf)
}


CalcForward <- function(y, A, B, istart) {
    # y             emitted symbols from an HMM
    # A, B, istart  HMM model to get likelihood of
    #
    # Returns:
    # L             log likelihood of the model
    # aft           forward probs for all timesteps, all states
    state.space <- 1:nrow(A)
    asf <- rep(0, length(y))
    aft <- matrix(rep(0, length(y)*length(state.space)), nrow=length(y))
    colnames(aft) <- state.space

    pf <- istart * B[, y[1]]
    sf <- sum(pf)
    pf <- pf / sf
    asf[1] <- sf
    aft[1, ] <- pf

    for (i in 2:length(y)) {
        nf <- sapply(state.space, function(x) {UpdateForward(pf, A, B, x, y[i])})
        sf <- sum(nf)
        asf[i] <- sf
        nf <- nf / sf
        aft[i, ] <- nf
        pf <- nf
    }

    L <- sum(log(asf))

    return(list(aft=aft, L=L))
}

UpdateBackward <- function(b, A, B, s, yi) {
    # calculate b_t-1 from b_t which is b argument
    nb <- sum(B[, yi] * A[s, ] * b)

    return(nb)
}

CalcBackward <- function(y, A, B) {
    # y  observations
    # A, B an HMM
    #
    # Returns:
    # abt  normalised backwards variables for all states all timesteps

    state.space <- 1:nrow(A)
    abt <- matrix(rep(0, length(y)*length(state.space)), nrow=length(y))
    colnames(abt) <- as.character(state.space)

    pb <- rep(1, length(state.space))
    pb <- pb / sum(pb)

    abt[length(y), ] <- pb

    for (i in (length(y)-1):1) {
        nb <- sapply(state.space, function(x) {UpdateBackward(pb, A, B, x, y[i+1])})
        nb <- nb / sum(nb)
        abt[i, ] <- nb
        pb <- nb
    }

    return(abt)
}

TestLikelihood <- function() {
    obs <- as.matrix(read.table("data/dsCG.dat"))
    obs <- as.vector(obs)

    trans.f <- "data/trans_hmm.dat"
    init.f <- "data/start_hmm.dat"
    emit.f <- "data/emit.dat"
    model <- ReadModel(trans.f, init.f, emit.f)

    f.res <- CalcForward(obs, model$A, model$B, model$istart)
    b.res <- CalcBackward(obs, model$A, model$B)

    return(f.res)
}

UpdateTransition <- function(A, B, y, aft, abt, s1, s2, L) {
    aft1 <- aft[1:(nrow(aft)-1), ]
    abt1 <- abt[2:(nrow(abt)), ]
    y1 <- y[2:length(y)]

    nt <- sum(aft1*abt1*A[s1, s2]*B[s2, y1])

    return(nt/L)
}

UpdateTransitionMatrix <- function(A, B, aft, abt, y, L) {
    # Returns:
    # nA new transition matrix

    nA <- matrix(rep(0, length(A)), nrow=nrow(model$A))
    state.space <- 1:nrow(A)

    for (i in state.space) {
        for (k in state.space) {
            nA[i, k] <- UpdateTransition(A, B, y, aft, abt, i, k, L)
        }
    }

    ni <- rowSums(nA)
    for (j in 1:nrow(nA)) {
        nA[j, ] <- nA[j, ] / ni[j]
    }

    return(nA)
}

UpdateEmission <- function(abt, aft, y, s, yi, L) {
    id <- rep(0, length(y))
    eq <- which(y == yi)
    neq <- setdiff(1:length(y), eq)
    id[eq] <- 1
    id[neq] <- 0

    ne <- sum(id * abt[, s] * aft[, s])

    return(ne/L)
}

UpdateEmissionMatrix <- function(B, aft, abt, y, L) {
    # Returns:
    # nB   new emission matrix
    nB <- matrix(rep(0, nrow(B)*ncol(B)), nrow=nrow(B))

    for (i in 1:nrow(nB)) {
        for (k in 1:ncol(nB)) {
            # i:state, k:symbol
            nB[i, k] <- UpdateEmission(abt, aft, y, i, k, L)
        }
    }

    #normalise, basically make probabilities out of counts
    ni <- rowSums(nB)
    for (j in 1:nrow(nB)) {
        nB[j, ] <- nB[j, ] / ni[j]
    }

    return(nB)
}

BaumWelch <- function(y, max.iter, A, B, istart) {
    # y observation
    # max.iter maximum number of iterations before stopping
    #
    # Returns:
    # A, B, istart  updated parameters

    # start by picking arbitrary parameters A, B, istart
    # for now take the A, B and istart given in the sheet
    f.res <- CalcForward(y, A, B, istart)
    print(f.res$L)
    for (i in 1:max.iter) {
        # at each iteration calculate aft, abt from A, B, istart
        # update A, B with UpdateEmissionMatrix, UpdateTransitionMatrix functions
        # repeat
        f.res <- CalcForward(y, A, B, istart)
        abt <- CalcBackward(y, A, B)
        aft <- f.res$aft

        A <- UpdateTransitionMatrix(A, B, aft, abt, y, f.res$L)
        B <- UpdateEmissionMatrix(B, aft, abt, y, f.res$L)

        f.res <- CalcForward(y, A, B, istart)
        print(f.res$L)
    }

    return(list(A=A, B=B, istart=istart))
}

UpdateViterbi <- function(v, A, B, yi, s) {
    nv <- B[s, yi] + max(v + A[, s])

    return(nv)
}

FindMaxState <- function(v, A, s) {
    i <- which.max(v + A[, s])

    return(i)
}

Viterbi <- function(y, A, B, istart) {
    A <- log(A)
    B <- log(B)
    state.space <- 1:nrow(A)
    ani <- matrix(rep(0, length(state.space)*length(y)), nrow=length(state.space))
    pv <- rep(1, length(state.space))

    for (i in 1:length(y)) {
        nv <- sapply(state.space, function(x) {UpdateViterbi(pv, A, B, y[i], x)})
        ni <- sapply(state.space, function(x) {FindMaxState(pv, A, x)})
        ani[, i] <- ni
        pv <- nv
    }


    ki <- which.max(nv)
    hs <- rep(0, length(y))
    hs[length(y)] <- ki

    for (j in (length(y)-1):1) {
        ki <- ani[ki, j+1]
        hs[j] <- ki
    }

    return(hs)
}


