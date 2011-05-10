bcsft.step <- function(rf = 0.5, s = 1, t = 0, type)
{
  if(missing(type)) {
    if(t == 0)
      type <- "bc"
    else if(s == 0)
      type <- "f2"
    else
      type <- "bcsft"
  }
  switch(type,
         f2 = mysteps.f2(rf, s, t),
         bc = mysteps.bc(rf, s, t),
         mysteps.default(rf, s, t))
}
mystep <- function(gen1 = 1, gen2 = 1, rf = 0.5, s = 1, t = 0)
{
  cross.scheme <- c(s,t)
  ret <- double(5)
  transpr <- double(10)
  z <- .C("step_wrap", as.integer(gen1), as.integer(gen2), as.double(rf), as.integer(cross.scheme),
          ret = as.double(ret), transpr = as.double(transpr))
  list(z$ret,z$transpr)
}
mysteps.default <- function(rf, s, t)
{
  out <- matrix(0,16,4)
  cross.type <- paste("BC", s, "F", t, sep = "")
  dimnames(out)[[2]] <- c("gen1","gen2", paste(cross.type, "b", sep = ""), cross.type)
  for(gen2 in 1:4) for(gen1 in 1:4) {
    tmp <- mystep(gen1,gen2,rf,s,t)
    if(gen1 == 1 & gen2 == 1)
      transpr <- tmp[[2]]
    out[gen2 + 4 * (gen1 - 1),] <- c(gen1, gen2, tmp[[1]][1:2])
  }
  dimnames(out)[[1]] <- paste(out[,1], out[,2], sep = ".")
  out[,-(1:2)] <- exp(out[,-(1:2)])

  transpr[8:10] <- exp(transpr[8:10])
  names(transpr) <- c("A1","B1","C","D","E","A0","B0","pr(1)","pr(2)","pr(3)")

  pr <- out
  mar <- c(transpr[8:10], NA)
  pr[,4] <- pr[,4] * mar[pr[,1]]
  mar <- mar[c(1:2,2:3)]
  pr[,3] <- pr[,3] * mar[pr[,1]]

  if(t == 0) {
    out <- out[out[,1] < 3 & out[,2] < 3,]
    pr <- pr[pr[,1] < 3 & pr[,2] < 3,]
    transpr <- transpr[c(1,2,4,8,9)]
  }
  else
    out[out[,1] > 3 | out[,2] > 3, 4] <- NA

  list(step = out, prob = pr, transpr = transpr)
}
mysteps.f2 <- function(rf, s, t)
{
  out <- matrix(0,16,6)
  cross.type <- paste("BC", s, "F", t, sep = "")
  dimnames(out)[[2]] <- c("gen1","gen2",paste(cross.type, "b", sep = ""),"F2b",cross.type,"F2")
  for(gen2 in 1:4) for(gen1 in 1:4) {
    tmp <- mystep(gen1,gen2,rf,s,t)
    if(gen1 == 1 & gen2 == 1)
      transpr <- tmp[[2]]
    out[gen2 + 4 * (gen1 - 1),] <- c(gen1, gen2, tmp[[1]][c(1,4,2,5)])
  }
  dimnames(out)[[1]] <- paste(out[,1], out[,2], sep = ".")
  out[,-(1:2)] <- exp(out[,-(1:2)])
  out[out[,1] > 3 | out[,2] > 3, 5:6] <- NA

  transpr[8:10] <- exp(transpr[8:10])
  names(transpr) <- c("A1","B1","C","D","E","A0","B0","pr(1)","pr(2)","pr(3)")

  list(step = out, transpr = transpr)
}
mysteps.bc <- function(rf, s, t)
{
  out <- matrix(0,4,5)
  cross.type <- paste("BC", s, "F", t, sep = "")
  dimnames(out)[[2]] <- c("gen1","gen2",paste(cross.type, "b", sep = ""),cross.type,"BC")
  for(gen2 in 1:2) for(gen1 in 1:2) {
    tmp <- mystep(gen1,gen2,rf,s,t)
    if(gen1 == 1 & gen2 == 1)
      transpr <- tmp[[2]]
    out[gen2 + 2 * (gen1 - 1),] <- c(gen1, gen2, tmp[[1]][1:3])
  }
  dimnames(out)[[1]] <- paste(out[,1], out[,2], sep = ".")
  out[,-(1:2)] <- exp(out[,-(1:2)])

  transpr[8:10] <- exp(transpr[8:10])
  names(transpr) <- c("A1","B1","C","D","E","A0","B0","pr(1)","pr(2)","pr(3)")
  transpr <- transpr[c(1,2,4,8,9)]

  list(step = out, transpr = transpr)
}

bcsft.prob <- function(rf, s, t, phase = FALSE)
{
  n.gen <- 2 + (t > 0) * (1 + phase)
  pr <- matrix(0, n.gen, n.gen)
  for(gen2 in seq(n.gen)) for(gen1 in seq(n.gen)) {
    tmp <- mystep(gen1, gen2, rf, s, t)
    if(gen1 == 1 & gen2 == 1)
      transpr <- tmp[[2]]
    pr[gen1,gen2] <- tmp[[1]][2 - phase]
  }
  ## Conditional prob(gen2|gen1).
  pr <- exp(pr)

  ## Joint prob(gen1,gen2).
  mar <- exp(transpr[8:10])
  if(phase)
    mar <- mar[c(1:2,2:3)]
  pr * mar[row(pr)]
}
