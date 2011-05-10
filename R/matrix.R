## genprobs: generate transition probabilites over generations.
## sumcts: summaries of expected counts over generations.
##
## both use transmat for transition matrices of bc and ft.
##
## F2 calcs.
## r 0.5 0.5
## t 2   3
## A 0   2/3
## B 1   4/3
## C 2   14/9
## D 0   2
## E 2   2

######################################################################################
## Transition probabilities.
transmat <- function(rf = 0.5)
{
  out <- list()
  out$ft <- transmat.ft(rf)
  out$bcs <- transmat.bcs(rf)
  out
}
transmat.ft <- function(rf = 0.5)
{
  ## Transition probability matrix for selfing.
  mat <- matrix(0, 10, 10)
  for(i in c(1,3,8,10))
    mat[i,i] <- 1
  tmp1 <- c(2,4,7,9)
  for(i in tmp1)
    mat[i,i] <- 0.5
  tmp2 <- c(1,3,4,1)
  tmp3 <- c(0,1,-1,0)
  for(i in seq(tmp1)) {
    mat[tmp1[i],tmp1[i] - tmp2[i]] <- 0.25
    mat[tmp1[i],tmp1[i] + tmp2[i] + tmp3[i]] <- 0.25
  }
  w <- 1 - rf
  mat[5,c(1,10)] <- mat[6,c(3,8)] <- w  ^ 2 / 4
  mat[6,c(1,10)] <- mat[5,c(3,8)] <- rf ^ 2 / 4
  mat[c(5,6),c(2,4,7,9)] <- rf * w / 2
  mat[5,5] <- mat[6,6] <- w  ^ 2 /2
  mat[5,6] <- mat[6,5] <- rf ^ 2 /2

  tmp1 <- c("A:AB.AB","B:AB.Ab","C:Ab.Ab","B:AB.aB","D:AB.ab","E:Ab.aB","B:Ab.ab","C:aB.aB","B:aB.ab","A:ab.ab")
  dimnames(mat) <- list(tmp1,tmp1)
  mat
}
transmat.bcs <- function(rf = 0.5)
{
  ## Transition probability matrix for selfing.
  mat <- matrix(0, 4, 4)
  ## Absorbing state: ab.ab.
  for(i in 4)
    mat[i,i] <- 1

  for(i in 2:3) {
    mat[i,i] <- 0.5
    mat[i,4] <- 0.5
  }
  
  w <- 1 - rf
  mat[1,c(1,4)] <- w / 2
  mat[1,2:3] <- rf / 2

  tmp1 <- c("AB.ab","Ab.ab","aB.ab","ab.ab")
  dimnames(mat) <- list(tmp1,tmp1)
  mat
}
genprob <- function(rf = 0.5, t = 2, init = pr,
                    mat = transmat(rf),
                    cross.scheme = c(0,t))
{
  ## Start with BCs
  s <- cross.scheme[1]
  if(s > 0) { ## BCs
    pr <- rep(0,4)
    pr[1] <- 1
    
    pr <- init
    
    for(i in seq(s))
      pr <- pr %*% mat$bcs
    
    ## For BCsFt, the initial state will depend on running BCs first.
    ## Translate from 4 BCs states to 10 Ft states.
    init <- rep(0,10)
    init[c(5,7,9,10)] <- pr
  }
  else { ## No BCs first.
    ## F1: all in state D:AB.ab
    pr <- rep(0,10)
    pr[5] <- 1
    ## This is the initial F1, so don't need to do it again.
    cross.scheme[2] <- cross.scheme[2] - 1
  }

  pr <- init
  t <- cross.scheme[2] ## Number of generations of selfing.
  if(t > 0) for(i in seq(t))
    pr <- pr %*% mat$ft

  pr <- c(pr)
  names(pr) <- dimnames(mat$ft)[[1]]
  pr
}
gennames <- function(cross.scheme)
{
  bc.gen <- cross.scheme[1]
  gen <- cross.scheme[2]
  gen.names <- NULL
  if(bc.gen > 0)
    gen.names <- paste("BC", seq(bc.gen), sep = "")
  if(gen > 0)
    gen.names <- c(gen.names, paste("BC", bc.gen, "F", seq(gen), sep = ""))
  gen.names
}  
genprobs <- function(rf = 0.5, gen = 10, bc.gen = 0, cross.scheme = c(bc.gen, gen))
{
  bc.gen <- cross.scheme[1]
  gen <- cross.scheme[2]
  gen.names <- gennames(cross.scheme)
  
  out <- matrix(NA, gen + bc.gen, 10)
  dimnames(out)[[1]] <- gen.names
  
  ## Transition matrices.
  mat <- transmat(rf)
  
  ## BCs.
  if(bc.gen > 0) {
    for(s in seq(cross.scheme[1])) {
      tmp <- genprob(rf, mat = mat, cross.scheme = c(s, 0))
      out[s,] <- tmp
    }
    tmp <- genprob(rf, 2, tmp, mat)
  }
  else
    tmp <- genprob(rf, 1, , mat)

  if(gen > 0) {
    out[bc.gen + 1, ] <- tmp
    if(gen > 1) for(i in seq(2, gen))
      out[bc.gen + i, ] <- tmp <- genprob(rf, 2, tmp, mat)
  }
  dimnames(out)[[2]] <- names(tmp)

  out
}
######################################################################################
## Counts
transct <- function()
{
  out <- list()
  out$ft <- transct.ft()
  out$bcs <- transct.bcs()
  out
}
transct.bcs <- function()
{
  ## Recombination count matrix for selfing.
  mat <- matrix(0, 4, 4)
  ## Absorbing state: ab.ab.
  for(i in 2:3)
    mat[1,i] <- 1

  tmp1 <- c("D:AB.ab","B:Ab.ab","B:aB.ab","A:ab.ab")
  dimnames(mat) <- list(tmp1,tmp1)
  mat
}
transct.ft <- function()
{
  ## Recombination count matrix for selfing.
  mat <- matrix(0, 10, 10)
  mat[6,c(1,10)] <- mat[5,c(3,8)] <- 2
  mat[c(5,6),c(2,4,7,9)] <- 1
  mat[5,6] <- mat[6,5] <- 2

  tmp1 <- c("A:AB.AB","B:AB.Ab","C:Ab.Ab","B:AB.aB","D:AB.ab","E:Ab.aB","B:Ab.ab","C:aB.aB","B:aB.ab","A:ab.ab")
  dimnames(mat) <- list(tmp1,tmp1)
  mat
}
multct <- function(t, init, mat, ct, numerator = rep(0, length(init)), denominator = init)
{
  if(t > 0) {
    matct <- mat * ct
  
    for(i in seq(t)) {
      ## Numerator takes previous numerator to new positions and adds new counts via denom %*% (mat*ct).
      numerator <- (numerator %*% mat) + (denominator %*% matct)
      denominator <- denominator %*% mat
    }
  }

  list(numerator = numerator, denominator = denominator)
}
sumct <- function(rf = 0.5, t = 2, init = pr, mat = transmat(rf), ct = transct(),
                  cross.scheme = c(0,t), numerator = FALSE)
{
  ## Start with BCs
  s <- cross.scheme[1]
  if(s > 0) { ## BCs
    pr <- rep(0,4)
    pr[1] <- 1
    
    pr <- init
    
    out <- multct(s, init, mat$bcs, ct$bcs)
    
    ## For BCsFt, the initial state will depend on running BCs first.
    ## Translate from 4 BCs states to 10 Ft states.
    init <- rep(0,10)
    init[c(5,7,9,10)] <- out$numerator
    out$numerator <- init
    init[c(5,7,9,10)] <- out$denominator
    out$denominator <- init
  }
  else { ## No BCs first.
    ## F1: all in state D:AB.ab
    pr <- rep(0,10)
    pr[5] <- 1

    pr <- init

    ## This is the initial F1, so don't need to do it again.
    cross.scheme[2] <- cross.scheme[2] - 1
    out <- multct(0, init)
  }

  t <- cross.scheme[2] ## Number of generations of selfing.
  out <- multct(t, init, mat$ft, ct$ft, out$numerator, out$denominator)
  if(numerator)
      out$numerator
  else {
    out <- out$numerator / out$denominator
    out[is.na(out)] <- 0
    out
  }
}
sumcts <- function(rf = 0.5, gen = 10, cross.scheme = c(0,gen), numerator = FALSE)
{
  bc.gen <- cross.scheme[1]
  gen <- cross.scheme[2]
  gen.names <- gennames(cross.scheme)
  
  out <- matrix(NA, gen + bc.gen, 10)
  dimnames(out)[[1]] <- gen.names

  ## Transition matrices.
  mat <- transmat(rf)
  ct <- transct()
  
  ## BCs.
  if(bc.gen > 0) {
    for(s in seq(cross.scheme[1])) {
      tmp <- sumct(rf, mat = mat, ct = ct, cross.scheme = c(s, 0), numerator = numerator)
      out[s,] <- tmp
    }
  }

  if(gen > 0) {
    for(i in seq(gen))
      out[bc.gen + i, ] <- tmp <- sumct(rf, mat = mat, ct = ct, cross.scheme = c(bc.gen, i),
                                        numerator = numerator)
    tmp <- dimnames(tmp)[[2]]
  }
  dimnames(out)[[2]] <- tmp
  out
}

######################################################################################
plotgen <- function(rf = 0.5, gen = 10, 
                    x = genprobs(rf, gen, cross.scheme = cross.scheme),
                    cross.scheme = c(0,gen), ylab = "prob", ...)
{
  tmp <- dimnames(x)[[2]]
  gen.names <- dimnames(x)[[1]]
  n.gen <- nrow(x)
  tmp1 <- substring(tmp, 1, 1)
  dat <- data.frame(state = ordered(rep(tmp, n.gen), tmp),
                    type = ordered(rep(tmp1, n.gen), sort(unique(tmp1))),
                    gen = ordered(rep(gen.names, rep(10, n.gen)), gen.names),
                    prob = c(t(x)))

  require(lattice)
  bc.gen <- cross.scheme[1]
  if(bc.gen > 0)
    panel.fun <- function(...) {
                    panel.stripplot(...)
                    panel.abline(v = bc.gen + 0.5, lty = 2, lwd = 2, col = "gray")
                  }
  else
    panel.fun <- panel.stripplot
  
  trellis.par.set(superpose.symbol = list(pch = 1:10))
  print(stripplot(prob ~ gen, dat, groups = state, type = "b", jitter = TRUE, ylab = ylab,
                  auto.key = list(space = "right"), scales=list(x=list(rot=45)),
                  panel = panel.fun))
  invisible(dat)
}
plotct <- function(rf = 0.5, gen = 10, 
                   x = sumcts(rf, gen, cross.scheme = cross.scheme, numerator = numerator),
                   cross.scheme = c(0,gen), ylab = "counts", ..., numerator = FALSE)
  plotgen(rf, gen, x, cross.scheme, ylab, ...)
