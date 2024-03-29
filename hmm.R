## ---- HMM ----
HMMGen <- function(A, emit, istart, N) {
    # A  transition matrix
    # emit  a function that given a hidden state outputs a symbol
    # istart  initial distribution
    # N  number of elements to output

    out.hs <- MarkovChain(A, istart, N)
    out.es <- sapply(out.hs, emit)

    return(list(hs=out.hs, es=out.es))
}

Emit <- function(B) {
    return( function(s) {
        sample(1:ncol(B), size=1, prob=B[s, ])
    })
}

PlotOut <- function(out) {
    # Expects out as given by HMMGen

    T <- 1:length(out$hs)
    states <- as.matrix(cbind(T, out$hs, out$es))

    ltyp <- rep("solid", ncol(states)-1)
    matplot(states[, 1], states[, 2:ncol(states)], type='l', ylab='State index',
            xlab = 'time', lwd=2, lty = ltyp, ylim = c(1, 6))

    legend("topright", c("Hidden States", "Emitted States"),
           col = 1:(ncol(states)-1), lty=ltyp)
}

## ---- ForwardAlg ----

CalcForward <- function(y, A, B, istart) {
    # y             observations
    # A, B, istart  HMM
    #
    # Returns:
    # L             log likelihood of the model defined by [A, B, istart]
    # aft           forward variables for all states, all timepoints
    # asf           scaling factors for all timepoints

    state.space <- 1:nrow(A)
    asf <- rep(0, length(y))
    aft <- matrix(rep(0, length(y)*length(state.space)), nrow=length(y))
    colnames(aft) <- state.space

    f <- istart * B[, y[1]]
    sf <- sum(f)
    f <- f / sf
    asf[1] <- sf
    aft[1, ] <- f
    for (i in 2:length(y)) {
        f <- sapply(state.space, function(s) { B[s, y[i]] * sum(f * A[, s])})
        sf <- sum(f)
        asf[i] <- sf
        f <- f / sf
        aft[i, ] <- f
    }

    L <- sum(log(asf))

    return(list(aft=aft, asf=asf, L=L))
}

# ---- BaumWelch ----

CalcBackward <- function(y, A, B, asf) {
    # y  observations
    # A, B HMM
    # asf  scaling factors from forward calculation
    #
    # Returns:
    # abt  normalised backwards variables for all states all timesteps

    state.space <- 1:nrow(A)
    abt <- matrix(rep(0, length(y)*length(state.space)), nrow=length(y))
    colnames(abt) <- as.character(state.space)

    b <- rep(1, length(state.space))
    abt[length(y), ] <- b
    for (i in (length(y)-1):1) {
        b <- sapply(state.space, function(s) {sum(B[, y[i+1]] * A[s, ] * b)})
        b <- b / asf[i+1]
        abt[i, ] <- b
    }

    return(abt)
}


UpdateParams <- function(A, B, istart, aft, abt, y) {
    n.states <- nrow(A)
    n.symb <- ncol(B)
    A1 <- matrix(rep(0, n.states**2), nrow=n.states)
    p.seq <- sum(aft[1, ] * abt[1, ])

    for (i in 1:n.states) {
        for (j in 1:n.states) {
            r <- sapply(1:(length(y)-2), function(t) {aft[t, i]*B[j, y[t+1]]*
                                                          A[i, j]*abt[t+2, j]})
            A1[i, j] <- sum(r) / p.seq
        }
    }

    A1 <- A1 / rowSums(A1)

    B1 <- matrix(rep(0, n.states*n.symb), nrow=n.states)
    for (i in 1:n.states) {
        gmi <- rep(0, length(y)-2)
        for (t in 1:(length(y)-2)){
            gm <- sum(sapply(1:n.states, function(j) {aft[t, i]*B[j, y[t+1]]*
                                                          A[i, j]*abt[t+2, j]}))
            gmi[t] <- gm
        }
        denom <- sum(gmi)
        for (k in 1:n.symb) {
            id <- rep(0, length(y)-2)
            eq <- which(y[1:(length(y)-2)] == k)
            id[eq] <- 1
            B1[i, k] <- sum(gmi * id) / denom
        }
    }

    t <- 1
    istart1 <- rep(0, n.states)
    for (i in 1:n.states) {
        init <- sum(sapply(1:n.states, function(j) {aft[t, i]*B[j, y[t+1]]*
                                                        A[i, j]*abt[t+2, j]}))
        istart1[i] <- init / p.seq
    }

    istart1 <- istart1 / sum(istart1)

    return(list(A=A1, B=B1, istart=istart1))
}

BaumWelch <- function(y, n.states, n.symb, max.iter, eps) {
    # Infer an HMM parameters from observations
    #
    # y             observations
    # n.states      number of states for the model
    # n.symb        number of symbols emitted
    # max.iter      maximum number of iterations
    # eps
    #
    # Returns:
    # [A, B, istart] inferred HMM from observations

    upt <- 1 / n.states
    upe <- 1 / n.symb
    s <- 0.2
    n.iter <- 0
    A <- matrix(abs(rnorm(n.states**2, upt, s)), nrow=n.states)
    A <- A / rowSums(A)

    B <- matrix(abs(rnorm(n.symb*n.states, upe, s)), nrow=n.states)
    B <- B / rowSums(B)

    istart <- rep(1/n.states, n.states)
    L0 <- 0
    L1 <- CalcForward(y, A, B, istart)$L
    dL <- abs(L1 - L0)
    while(n.iter < max.iter &&
          dL > eps) {
        L0 <- L1
        f.res <- CalcForward(y, A, B, istart)
        aft <- f.res$aft
        abt <- CalcBackward(y, A, B, f.res$asf)

        nParams <- UpdateParams(A, B, istart, aft, abt, y)
        A <- nParams$A
        B <- nParams$B
        istart <- nParams$istart
        L1 <- CalcForward(y, A, B, istart)$L
        dL <- abs(L1 - L0)
        n.iter <- n.iter + 1
    }

    return(list(A=A, B=B, istart=istart))
}

## ---- Viterbi -----
Viterbi <- function(y, A, B, istart) {
    # y  observations
    # A, B, istart  HMM
    A <- log(A)
    B <- log(B)
    state.space <- 1:nrow(A)
    ani <- matrix(rep(0, length(state.space)*length(y)), nrow=length(state.space))
    v <- rep(1, length(state.space))

    for (i in 1:length(y)) {
        v <- sapply(state.space, function(s) { B[s, y[i]] + max(v + A[, s]) })
        ni <- sapply(state.space, function(s) { which.max(v + A[, s]) })
        ani[, i] <- ni
    }

    ki <- which.max(v)
    hs <- rep(0, length(y))
    hs[length(y)] <- ki

    for (j in (length(y)-1):1) {
        ki <- ani[ki, j+1]
        hs[j] <- ki
    }

    return(hs)
}

## ---- Test ----
PlotContent <- function(vhs, content) {
    T <- 1:length(content)
    dat <- as.matrix(cbind(hs, cg))

    ltyp <- rep("solid", ncol(dat)-1)
    matplot(T, dat, type='l', ylab='State index/CG content',
            xlab = 'time', lwd=2, lty = ltyp, ylim=c(0, 3))

    legend("topright", c("Hidden States", "GC content(%)"),
           col = 1:(ncol(dat)), lty=ltyp)
}

TestLikelihood <- function() {
    obs <- as.matrix(read.table("data/dsCG.dat"))
    obs <- as.vector(obs)

    ptrans.f <- "data/trans_hmm.dat"
    init.f <- "data/start_hmm.dat"
    emit.f <- "data/emit.dat"
    model <- ReadModel(trans.f, init.f, emit.f)

    f.res <- CalcForward(obs, model$A, model$B, model$istart)
    b.res <- CalcBackward(obs, model$A, model$B)

    return(f.res)
}



