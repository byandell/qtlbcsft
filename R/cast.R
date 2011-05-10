cast <- function(rf = 0.5, gen = 10, trans = calctrans(rf, gen))
{
## Efficient calcs here to check.

  ncount <- matrix(0, gen, 8)
  dimnames(ncount) <- list(seq(gen), c("NDDA","NDDB","NDDC","NDDD","NDBA","NDBB","NDBA2","NDBB2"))
  count <- ncount

  ## PDx = probability of getting to x from D at generation t.
  ## PBA = probability of getting to A from B in t generations.
  ## MDx = probability * count of recombinations in t generations and stay in D or E.
  ## Nxyz = probability * count ending at z passing through x and later y.
  ## expected count is sum of Nxyz divided by probability for end state.

  tmp <- c("t","sbeta1","sgamma1","sbeta2","sgamma2", "s2beta1", "s2gamma1", "s2beta2", "s2gamma2","k1b","k1g","k2b","k2g","k1b2","k1g2","k2b2","k2g2")
  diag <- matrix(0, gen, length(tmp))
  dimnames(diag) <- list(seq(gen), tmp)
  for(t in seq(2,gen)) {
    t1 = t - 1.0;
    t2 = R_pow(2.0, -t1);

    r2 = rf * rf;
    w2 = (1.0 - rf) * (1.0 - rf);
    rw = rf * (1.0 - rf);

    alpha = r2 / (w2 + r2);
    beta = (w2 + r2) / 2.0; 
    gamma = (w2 - r2) / 2.0;

    beta1 = R_pow(beta, t1);
    gamma1 = 1.0
    gamma2 = 1.0
    if(t > 1) gamma1 = R_pow(gamma, t1);
    if(t > 2) gamma2 = R_pow(gamma, t1 - 1);

    sbeta1 = (1.0 - beta1) / (1.0 - beta); ## SFt
    sbeta2 = (1.0 - beta1 / beta) / (1.0 - beta); ## SF(t-1)
    s2beta1 = (t2 - beta1) / (1.0 - 2.0 * beta);
    s2beta2 = (2.0 * t2 - (beta1 / beta)) / (1.0 - 2.0 * beta);

    sgamma1 = 1
    sgamma2 = 0
    s2gamma1 = t2
    s2gamma2 = 0
    if(t > 1) {
      sgamma2 = 1
      s2gamma2 = t2 * 2
    }
    if(gamma > 0) {
      sgamma1 = (1.0 - gamma1) / (1.0 - gamma); ## SGt
      sgamma2 = (1.0 - gamma2) / (1.0 - gamma); ## SG(t-1)
      s2gamma1 = (t2  - gamma1) / (1.0 - 2.0 * gamma);
      s2gamma2 = (t2 * 2.0 - (gamma1 / gamma)) / (1.0 - 2.0 * gamma);
    }

    ## Need to get the *2s to line up properly.
    k1b <- kptothek(t-1, beta) / beta
    k2b <- t2 * kptothek(t-1, 2.0 * beta) / (2 * beta)
    k1g <- 0.0
    k2g <- 0.0
    k1g2 <- 0.0
    k2g2 <- 0.0
    k1b2 <- 0.0
    k2b2 <- 0.0
    if(t > 2) {
      k1g <- 1.0
      k2g <- t2 ## not sure
      if(t > 3) {
        k1g2 <- 1.0
        k2g2 <- 2 * t2
      }
      k1b2 <- kptothek(t-2, beta) / beta
      k2b2 <- 2 * t2 * kptothek(t-2, 2.0 * beta) / (2 * beta)
    }
    if(gamma > 0) {
      k1g <- kptothek(t-1, gamma) / gamma
      k2g <- t2 * kptothek(t-1, 2.0 * gamma) / (2.0 * gamma)
      k1g2 <- kptothek(t-2, gamma) / gamma
      k2g2 <- 2 * t2 * kptothek(t-2, 2.0 * gamma) / (2.0 * gamma)
    }

    diag[t,] <- c(t, sbeta1, sgamma1, sbeta2, sgamma2, s2beta1, s2gamma1, s2beta2, s2gamma2, k1b,k1g,k2b,k2g, k1b2,k1g2,k2b2,k2g2)
    
    count[t,"NDDA"] <- sum(trans[seq(t-1),"PDD"]) ## OK
    count[t,"NDDB"] <- sum(trans[seq(t-1),"PDD"] * trans[seq(t-1),"MDD"]) ## OK
    count[t,"NDDC"] <- sum(trans[seq(t-1),"PDD"] * 2 ^ -seq(t-2,0)) ## OK
    count[t,"NDDD"] <- sum(trans[seq(t-1),"PDD"] * 2 ^ -seq(t-2,0) * trans[seq(t-1),"MDD"])

    if(t > 2) {
      count[t,"NDBA"] <- sum(trans[seq(t-2),"PDD"] * trans[seq(t-2,1),"PBA"])
      count[t,"NDBB"] <- sum(trans[seq(t-2),"PDD"] * trans[seq(t-2,1),"PBA"] * trans[seq(t-2),"MDD"])
    }

    ncount[t,"NDDA"] <- 0.5 * (sbeta1 + sgamma1)
    ncount[t,"NDDB"] <- r2 * 0.5 * (k1b - k1g)
    ncount[t,"NDDC"] <- (s2beta1 + s2gamma1) ## wrong
    ncount[t,"NDDD"] <- 2 * r2 * (k2b - k2g)
    if(t > 1) {
      count[t,"NDBA2"] <- 0.5 * (count[t-1,"NDDA"] - count[t-1,"NDDC"])
      count[t,"NDBB2"] <- 0.5 * (count[t-1,"NDDB"] - count[t-1,"NDDD"])
      ncount[t,"NDBA2"] <- 0.5 * (ncount[t-1,"NDDA"] - ncount[t-1,"NDDC"])
      ncount[t,"NDBB2"] <- 0.5 * (ncount[t-1,"NDDB"] - ncount[t-1,"NDDD"])
      ncount[t,"NDBA"] <- 0.25 * ((sbeta2 + sgamma2) - (s2beta2 + s2gamma2)) ## wrong
      ncount[t,"NDBB"] <- 0.5 * r2 * (0.5 * (k1b2 - k1g2) - (k2b2 - k2g2))
    }
  }

  lapply(list(diag = diag, count = count, ncount = ncount), round, 5)
}
ftcast <- function(rf = 0.5, gen = 10, trans = calctrans(rf, gen))
{
## Efficient calcs here to check.

  ocounts <- calccounts(rf, gen, trans)
  prob <- calcprobs(rf, gen, trans, TRUE)

  count <- matrix(0, gen, 12)
  dimnames(count) <- list(seq(gen),
                          c("NDDA","NDEA","NDBA","NEBA","NBA",
                            "NDDB","NDEB",
                            "NDDC","NDEC","NDBC","NEBC","NBC"))
  counts <- matrix(0, gen, 5)
  dimnames(counts) <- list(seq(gen), LETTERS[1:5])

  ## PDx = probability of getting to x from D at generation t.
  ## PBA = probability of getting to A from B in t generations.
  ## MDx = probability * count of recombinations in t generations and stay in D or E.
  ## Nxyz = probability * count ending at z passing through x and later y.
  ## expected count is sum of Nxyz divided by probability for end state.

  diag <- matrix(0, gen, 11)
  dimnames(diag) <- list(seq(gen),
                         c("t","beta","gamma","beta1","gamma1","sbeta1","sgamma1","s2beta1","s2gamma1","k2b","k2g"))
  for(t in seq(2,gen)) {
    t1 = t - 1.0;
    t2 = R_pow(2.0, -t1);

    r2 = rf * rf;
    w2 = (1.0 - rf) * (1.0 - rf);
    rw = rf * (1.0 - rf);

    alpha = r2 / (w2 + r2);
    beta = (w2 + r2) / 2.0; 
    gamma = (w2 - r2) / 2.0;

    beta1 = R_pow(beta, t1);
    gamma1 = 1.0
    gamma2 = 1.0
    if(t > 1) gamma1 = R_pow(gamma, t1);
    if(t > 2) gamma2 = R_pow(gamma, t1 - 1);

    sbeta1 = (1.0 - beta1) / (1.0 - beta); ## SFt
    sbeta2 = (1.0 - beta1 / beta) / (1.0 - beta); ## SF(t-1)
    s2beta1 = (t2 - beta1) / (1.0 - 2.0 * beta);
    s2beta2 = (2 * t2 - (beta1 / beta)) / (1.0 - 2.0 * beta);

    sgamma1 = 0
    sgamma2 = 0
    s2gamma1 = 0
    s2gamma2 = 0
    if(gamma > 0) {
      sgamma1 = (gamma - gamma1) / (1.0 - gamma); ## SGt
      sgamma2 = (gamma - gamma2) / (1.0 - gamma); ## SG(t-1)
      s2gamma1 = (t2 * 2.0 * gamma - gamma1) / (1.0 - 2.0 * gamma);
      s2gamma2 = (t2 * 4.0 * gamma - (gamma1 / gamma)) / (1.0 - 2.0 * gamma);
    }

    k1b <- kptothek(t-1, beta) / beta
    k2b <- t2 * kptothek(t-1, 2.0 * beta) / (2 * beta)
    k1g <- 0.0
    k2g <- 0.0
    k1g2 <- 0.0
    k2g2 <- 0.0
    k1b2 <- 0.0
    k2b2 <- 0.0
    if(t > 2) {
      k1g <- 1.0
      k2g <- t2
       if(t > 3) {
        k1g2 <- 1.0
        k2g2 <- 2 * t2
      }
     k1b2 <- kptothek(t-2, beta) / beta
      k2b2 <- 2 * t2 * kptothek(t-2, 2.0 * beta) / (2 * beta)
    }
    if(gamma > 0) {
      k1g <- kptothek(t-1, gamma) / gamma
      k2g <- t2 * kptothek(t-1, 2.0 * gamma) / (2.0 * gamma)
      k1g2 <- kptothek(t-2, gamma) / gamma
      k2g2 <- 2 * t2 * kptothek(t-2, 2.0 * gamma) / (2.0 * gamma)
    }

    diag[t,] <- c(t, beta, gamma, beta1, gamma1, sbeta1, sgamma1, s2beta1, s2gamma1,k2b,k2g)
    
    count[t,"NDDA"] <- (w2 / 4) * (r2 / 2) * (k1b - k1g) ## OK
    count[t,"NDDB"] <- (rw / 2) * (r2 * (k2b - k2g) + 0.5 * (s2beta1 + s2gamma1)) ## OK
    count[t,"NDDC"] <- (r2 / 4) * ((r2 / 2) * (k1b - k1g) + (sbeta1 + sgamma1)) ## OK
    
    if(t > 2) {
      count[t,"NDEA"] <- (r2 / 4) * ((r2 / 2) * (k1b + k1g) + (sbeta1 - sgamma1)) ## OK
      count[t,"NDEB"] <- (rw / 2) * (r2 * (k2b + k2g) + 0.5 * (s2beta1 - s2gamma1)) ## OK
      count[t,"NDEC"] <- (w2 / 4) * (r2 / 2) * (k1b + k1g) ## OK

      count[t,"NDBA"] <- rw * (0.25 * ((sbeta2 + sgamma2) - (s2beta2 + s2gamma2)) +
                               0.5 * r2 * (0.5 * (k1b2 - k1g2) - (k2b2 - k2g2)))
      count[t,"NDBC"] <- rw * (0.25 * ((sbeta2 + sgamma2) - (s2beta2 + s2gamma2)) +
                               0.5 * r2 * (0.5 * (k1b2 - k1g2) - (k2b2 - k2g2)))
      if(t > 3) {
        count[t,"NEBA"] <- rw * (0.25 * ((sbeta2 - sgamma2) - (s2beta2 - s2gamma2)) +
                               0.5 * r2 * (0.5 * (k1b2 + k1g2) - (k2b2 + k2g2)))
        count[t,"NEBC"] <- rw * (0.25 * ((sbeta2 - sgamma2) - (s2beta2 - s2gamma2)) +
                               0.5 * r2 * (0.5 * (k1b2 + k1g2) - (k2b2 + k2g2)))
      }
    }
  }

  count[, "NBA"] <- apply(count[,c("NDBA","NEBA")], 1, sum)
  count[, "NBC"] <- apply(count[,c("NDBC","NEBC")], 1, sum)

  counts[, "A"] <- apply(count[,c("NDDA","NDEA","NDBA","NEBA")], 1, sum)
  counts[, "B"] <- apply(count[,c("NDDB","NDEB")], 1, sum)
  counts[, "C"] <- apply(count[,c("NDDC","NDEC","NDBC","NEBC")], 1, sum)

##  sapply(list(diag = diag, count = cbind(count[,c("NDDA","NDDC","NDEA","NDEC")],
##        ocount[,c("NDDA","NDDC","NDEA","NDEC")])), round, 5)
##  tmp <- count[,c("NDBA","NDBC","NEBA","NEBC")]
##  tmp1 <- ocount[,c("NDBA","NDBC","NEBA","NEBC")]
##  sapply(list(diag = diag, count = cbind(tmp, NDB=apply(tmp,1,sum),
##                             tmp1, ODB=apply(tmp1,1,sum))), round, 5)

  list(counts = counts, ocounts = ocounts)
}
