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
    trans.f <- "data/trans.dat"
    init.f <- "data/start_hmm1.dat"
    emit.f <- "data/emit1.dat"
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

PlotContent <- function(vhs, content) {
    T <- 1:length(content)
    dat <- as.matrix(cbind(hs, cg))

    ltyp <- rep("solid", ncol(dat)-1)
    matplot(T, dat, type='l', ylab='State index',
            xlab = 'time', lwd=2, lty = ltyp)

    legend("topright", c("Hidden States", "GC content(%)"),
           col = 1:(ncol(dat)-1), lty=ltyp)
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

Ksi <- function(m, i, j, A, B, aft, abt, y) {
    # Calculate ksi var for transition i->j
    # at time m
    nf <- sum(aft[m, ] * abt[m, ])
    k <- (aft[m, i]*abt[m+2, j]*B[j, y[m+1]]*A[i, j]) / nf

    return(k)
}

GammaH <- function(m, i, A, B, aft, abt, y) {
    n.states <- nrow(A)
    return(sum(sapply(1:n.states, function(x) {Ksi(m, i, x, A, B, aft, abt, y)})))
}

UpdateTransition <- function(A, B, aft, abt, y) {
    n.states <- nrow(A)
    n.obs <- length(y)
    nA <- matrix(rep(0, n.states*n.states), nrow=n.states)
    G <- rep(0, n.states)

    for (i in 1:n.states) {
        gammaij <- sum(sapply(1:(n.obs-2), function(m) {GammaH(m, i, A, B, aft,
                                                               abt, y)}))
        G[i] <- gammaij
        for (j in 1:n.states) {
            ksiij <- sum(sapply(1:(length(y)-2), function(m) {Ksi(m, i, j, A,
                                                                  B, aft,
                                                                  abt, y)}))
            nA[i, j] <- ksiij / gammaij
        }
    }

    return(list(nA=nA, G=G))
}

UpdateStart <- function(A, B, aft, abt, y) {
    n.states <- nrow(A)
    start <- sapply(1:n.states, function(x) {GammaH(1, x, A, B, aft, abt, y)})
    start <- start / sum(start)

    return(start)
}

UpdateEmission <- function(A, B, aft, abt, y, G) {
    n.states <- nrow(A)
    n.obs <- length(y)
    n.symb <- ncol(B)

    nB <- matrix(rep(0, n.states*n.symb), nrow=n.states)

    for (j in 1:n.states) {
        gni <- sapply(1:(n.obs-2), function(x) {GammaH(x, j, A, B, aft, abt, y)})
        for (n in 1:n.symb) {
            id <- rep(0, n.obs-2)
            eq <- which(y[1:(n.obs-2)] == n)
            id[eq] <- 1

            nB[j, n] <- sum(gni * id)
        }
    }

    ni <- rowSums(nB)
    for (j in 1:nrow(nB)) {
        nB[j, ] <- nB[j, ] / ni[j]
    }

    return(nB)
}

PickInitEmission <- function(y, n.states, n.symb) {
    B <- matrix(rep(0, n.states*n.symb), nrow=n.states)

    kres <- kmeans(y, n.states)
    print(kres$centers)
    for (i in 1:n.states) {
        cl <- y[which(kres$cluster == i)]
        p <- table(cl) / sum(table(cl))
        B[i, as.integer(names(p))] <- p
    }

    return(B)
}

BaumWelch <- function(y, n.states, n.symb, max.iter, A, B, eps = 0.01) {
    n.iter <- 0
    L0 <- 0
    #pick arbitrary model parameters A, B, istart
    up <- 1 / n.states
    #A <- matrix(rep(up, n.states**2), nrow=n.states)
    #A <- matrix(c(0.8, 0.2, 0.1, 0.9), nrow=2, byrow = T)
    #B <- PickInitEmission(y, n.states, n.symb)
    istart <- rep(up, n.states)
    L1 <- CalcForward(y, A, B, istart)$L
    dL <- L1 - L0

    while(abs(dL) > eps &&
          n.iter < max.iter) {
        L0 <- L1
        aft <- CalcForward(y, A, B, istart)$aft
        abt <- CalcBackward(y, A, B)

        upd.res <- UpdateTransition(A, B, aft, abt, y)
        A <- upd.res$nA
        B <- UpdateEmission(A, B, aft, abt, y, upd.res$G)
        #istart <- UpdateStart(A, B, aft, abt, y)
        L1 <- CalcForward(y, A, B, istart)$L
        print(L1)
        dL <- L1 - L0
        n.iter <- n.iter + 1
    }

    return((list(A=A, B=B, istart=istart)))


    ## for (i in 1:50) {
    ##     f.res <- CalcForward(y, A, B, istart)
    ##     print(f.res$L)
    ##     abt <- CalcBackward(y, A, B)
    ##     aft <- f.res$aft

    ##     A <- UpdateTransition(list(A=A, B=B, istart=istart), aft, abt, y, 2)
    ##     B <- UpdateEmission(list(A=A, B=B, istart=istart), aft, abt, y)
    ##     istart <- UpdateStart(list(A=A, B=B, istart=istart), aft, abt, y, 2)
    ## }
    ## return(list(A=A, B=B, istart=istart))
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


