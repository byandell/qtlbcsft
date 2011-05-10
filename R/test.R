mytest <- function(obs1 = 1, obs2 = 1, rf = 0.5, s = 1, t = 0)
{
  cross.scheme <- c(s,t)
  init <- double(6)
  stepb <- nrec <- double(32)
  emit <- step <- double(18)
  transexp <- transpr <- double(10)
  z <- .C("bcsft_wrap", rf = as.double(rf), cross.scheme = as.integer(cross.scheme),
          init = as.double(init), emit = as.double(emit), step = as.double(step),
          stepb = as.double(stepb), nrec = as.double(nrec),
          transpr = as.double(transpr), transexp = as.double(transexp))
  z$init <- matrix(z$init, 3, 2)
  z$emit <- array(z$emit, c(3, 3, 2))
  z$step <- array(z$step, c(3, 3, 2))
  z$stepb <- array(z$stepb, c(4, 4, 2))
  z$nrec <- array(z$nrec, c(4, 4, 2))
  if(t == 0) {
    z$init <- z$init[1:2,]
    z$emit <- z$emit[1:2,1:2,]
    z$step <- z$step[1:2,1:2,]
    z$stepb <- z$stepb[1:2,1:2,]
    z$nrec <- z$nrec[1:2,1:2,]
  }
  z
}

######################################################################################

######################################################################################
bcsft.init <- function(s = 1, t = 0)
{
  cross.scheme <- c(s,t)
  out <- matrix(NA, 5, 4)
  gens <- 1:4
  if(t == 0)
    gens <- 1:2
  for(gen in gens) {
    z <- .C("init_wrap", as.integer(gen), as.integer(cross.scheme),
            ret = as.double(rep(0, 5)))
    out[,gen] <- z$ret
  }
  dimnames(out) <- list(c("BCsFtb","F2b","BCsFt","F2","BC"), 1:4)
  if(t == 0)
    out <- out[c(1,3,5),1:2]
  else if(s == 0)
    out <- out[-5,]
  out[out == 0] <- NA
  exp(out)
}

######################################################################################
bcsft.nrec <- function(rf = 0.5, s = 0, t = 2)
{
  cross.scheme <- c(s,t)
  out <- matrix(0,16,5)
  bcsft <- paste("BC",s,"F",t,sep = "")
  dimnames(out)[[2]] <- c("gen1","gen2",bcsft,"F2","BC")
  for(gen2 in 1:4) for(gen1 in 1:4) {
    ret <- double(3)
    z <- .C("nrec_wrap", as.integer(gen1), as.integer(gen2), as.double(rf), as.integer(cross.scheme),
            ret = as.double(ret))
    out[gen2 + 4 * (gen1 - 1),] <- c(gen1, gen2, z$ret[1:3])
  }
  dimnames(out)[[1]] <- paste(out[,1], out[,2], sep = ".")
  if(t == 0)
    out <- out[out[,1] < 3 & out[,2] < 3, ]
  attr(out, "scheme") <- cross.scheme
  
  out
}

######################################################################################
