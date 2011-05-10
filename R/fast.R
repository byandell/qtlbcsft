ftfast <- function(rf = 0.5, gen = 10, trans = calctrans(rf, gen), detail = FALSE)
{
## Efficient calcs here to check.

  ocount <- calccounts(rf, gen, trans, TRUE)
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
    beta1 = R_pow(beta, t1);
    beta2 = 1.0;
    if(t > 2.0) beta2 = beta1 / beta;
    sbeta1 = (1.0 - beta1) / (1.0 - beta); ## SFt
    sbeta2 = (1.0 - beta1 / beta) / (1.0 - beta); ## SF(t-1)
    s2beta1 = (t2 - beta1) / (1.0 - 2.0 * beta);
    s2beta2 = (2 * t2 - (beta1 / beta)) / (1.0 - 2.0 * beta);

    gamma = (w2 - r2) / 2.0;
    gamma1 = 1.0
    if(t > 1) gamma1 = R_pow(gamma, t1);
    gamma2 = 1.0
    if(t > 2) gamma2 = R_pow(gamma, t1 - 1);
    sgamma1 = 1
    sgamma2 = 0
    s2gamma1 = t2
    s2gamma2 = 0
    if(t > 1) {
      sgamma2 = 1
      s2gamma2 = t2 * 2
    }
    if(gamma > 0) {
      sgamma1 = (1.0- gamma1) / (1.0 - gamma); ## SGt
      sgamma2 = (1.0 - gamma2) / (1.0 - gamma); ## SG(t-1)
      s2gamma1 = (t2 - gamma1) / (1.0 - 2.0 * gamma);
      s2gamma2 = (2 * t2 - (gamma1 / gamma)) / (1.0 - 2.0 * gamma);
    }

    ## kptothek(t,p) = sum of k * p^k from 1 to t-1.
    k1b <- kptothek(t-1, beta, beta1) / beta
    k2b <- t2 * kptothek(t-1, 2.0 * beta, beta1 / t2) / (2 * beta)
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

    ## Numerator for expected counts
    ## End point is B. Already have this.
    count[t,"NDDB"] <- rw * (r2 * (k2b - k2g) + 0.5 * (s2beta1 + s2gamma1)) ## OK
    if(t > 2)
      count[t,"NDEB"] <- rw * (r2 * (k2b + k2g) + 0.5 * (s2beta1 - s2gamma1)) ## OK

    ## End point is A or C. Newly revised.
    ndda <- (r2 / 2) * (k1b - k1g)
    count[t,"NDDA"] <- (w2 / 4) * ndda ## OK
    count[t,"NDDC"] <- (r2 / 4) * (ndda + (sbeta1 + sgamma1)) ## OK

    
    if(t > 2) {
      ndea <- (r2 / 2) * (k1b + k1g)
      count[t,"NDEA"] <- (r2 / 4) * (ndea + (sbeta1 - sgamma1)) ## OK
      count[t,"NDEC"] <- (w2 / 4) * ndea ## OK

      nbab <- rw * (0.25 * (sbeta2 - s2beta2) + 0.5 * r2 * (0.5 * k1b2 - k2b2))
      nbag <- rw * (0.25 * (sgamma2 - s2gamma2) - 0.5 * r2 * (0.5 * k1g2 - k2g2))
      count[t,"NDBA"] <- nbab + nbag
      if(t > 3)
        count[t,"NEBA"] <- nbab - nbag
    }
  }

  counts[, "B"] <- apply(count[,c("NDDB","NDEB")], 1, sum)
  
  counts[, "A"] <- apply(count[,c("NDDA","NDEA","NDBA","NEBA")], 1, sum)
  counts[, "C"] <- apply(count[,c("NDDC","NDEC","NDBA","NEBA")], 1, sum)

  if(detail)
    count
  else
    counts
}
