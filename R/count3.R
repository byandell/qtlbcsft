## calcprobs and calccounts use calctrans.
## Algorithms similar to what was translated into C.

###############################################################################
ftcount <- function(rf = 0.5, gen = 10)
{
## Now have all the counts and probs correct.
## Need to put efficient calc here to check.

  ## expected count is sum of Nxyz divided by probability for end state.

  trans <- calctrans(rf, gen)

  probs <- calcprobs(rf, gen, trans)
  counts <- calccounts(rf, gen, trans)

  ## Expected counts.
  expct <- counts / probs
  expct[probs == 0] <- 0

  list(trans = trans, probs = probs, counts = counts, expct = expct)
  expct
}
###############################################################################
calctrans <- function(rf = 0.5, gen = 10)
{
  ## Transitions between states.
  ## PDx = probability of getting to x from D at generation t.
  ## PBA = probability of getting to A from B in t generations.
  ## MDx = probability * count of recombinations in t generations and stay in D or E.
  ## SDx = alternative to MDx.
  trans <- matrix(0, gen, 7)
  dimnames(trans) <- list(seq(gen), c("PDD","PDE","PBA","MDD","MDE","SDD","SDE"))

  for(t in seq(gen)) {
    t1 = t - 1.0;
    t2 = R_pow(2.0, -t1);

    r2 = rf * rf;
    w2 = (1.0 - rf) * (1.0 - rf);
    rw = rf * (1.0 - rf);

    ## Need better notation here. Want to eliminate alpha?
    alpha = r2 / (w2 + r2);
    beta = (w2 + r2) / 2.0; 
    gamma = (w2 - r2) / 2.0;
    beta1 = 1.0;
    gamma1 = 1.0;
    beta2 = beta1;
    gamma2 = 1.0;
    if(t > 1.0) {
      beta1 = R_pow(beta, t1);
      gamma1 = R_pow(gamma, t1);
    }
    if(t > 2.0) {
      beta2 = beta1 / beta;
      gamma2 = 0;
      if(gamma > 0) gamma2 = gamma1 / gamma;
    }

    alpha21 = 1.0 - 2.0 * alpha;
    alpha2 = 1.0;
    if(t > 2.0) alpha2 = R_pow(alpha21, t - 2.0);
    alpha1 = alpha2 * alpha21;

    ## now the terms we really want!
    trans[t,"PDD"] <- (beta1 + gamma1) / 2.0;
    trans[t,"PDE"] <- (beta1 - gamma1) / 2.0;
    trans[t,"PBA"] <- 2 * (1 - 2 ^ -t) / 4 ## 1 to t-1 steps in B, then A

    trans[t,"MDD"] <- 2.0 * t1 * alpha * (1.0 - alpha2) / (1.0 + alpha1); ## D->D
    trans[t,"MDE"] <- 2.0 * t1 * alpha * (1.0 + alpha2) / (1.0 - alpha1); ## E->E

    trans[t,"SDD"] <- t1 * r2 * (beta2 - gamma2) / (beta1 + gamma1); ## D->D
    trans[t,"SDE"] <- t1 * r2 * (beta2 + gamma2) / (beta1 - gamma1); ## E->E
  }
  trans
}
###############################################################################
calcprobs <- function(rf = 0.5, gen = 10, trans = calctrans(rf, gen),
                      detail = FALSE)
{
  ## Probabilities to end states.
  ## Nxyz = probability ending at z passing through x and later y.

  r2 = rf * rf;
  w2 = (1.0 - rf) * (1.0 - rf);
  rw = rf * (1.0 - rf);

  prob <- matrix(0, gen, 12)
  dimnames(prob) <- list(seq(gen),
                          c("NDDA","NDEA","NDBA","NEBA","NBA",
                            "NDDB","NDEB",
                            "NDDC","NDEC","NDBC","NEBC","NBC"))
  probs <- matrix(0, gen, 5)
  dimnames(probs) <- list(seq(gen), LETTERS[1:5])

  probs[, c("D","E")] <- trans[, c("PDD","PDE")]
  if(gen > 1) for(t in seq(2,gen)) {
    prob[t,"NDDA"] <- sum(trans[seq(t-1),"PDD"] * (w2 / 4))
    prob[t,"NDDB"] <- sum(trans[seq(t-1),"PDD"] * (rw / 2) * 2 ^ -seq(t-2,0))
    prob[t,"NDDC"] <- sum(trans[seq(t-1),"PDD"] * (r2 / 4))
    if(t > 2) {
      prob[t,"NDEA"] <- sum(trans[seq(t-1),"PDE"] * (r2 / 4))
      prob[t,"NDEB"] <- sum(trans[seq(t-1),"PDE"] * (rw / 2) * 2 ^ -seq(t-2,0))
      prob[t,"NDEC"] <- sum(trans[seq(t-1),"PDE"] * (w2 / 4))
      prob[t,"NDBA"] <- 2 * sum(trans[seq(t-2),"PDD"] * (rw / 2) * trans[seq(t-2,1),"PBA"])
      prob[t,"NDBC"] <- 2 * sum(trans[seq(t-2),"PDD"] * (rw / 2) * trans[seq(t-2,1),"PBA"])
      if(t > 3) {
        prob[t,"NEBA"] <- 2 * sum(trans[seq(t-2),"PDE"] * (rw / 2) * trans[seq(t-2,1),"PBA"])
        prob[t,"NEBC"] <- 2 * sum(trans[seq(t-2),"PDE"] * (rw / 2) * trans[seq(t-2,1),"PBA"])
      }
    }
  }
  prob[, "NBA"] <- apply(prob[,c("NDBA","NEBA"), drop = FALSE], 1, sum)
  prob[, "NBC"] <- apply(prob[,c("NDBC","NEBC"), drop = FALSE], 1, sum)
  probs[, "A"] <- apply(prob[,c("NDDA","NDEA","NDBA","NEBA"), drop = FALSE], 1, sum)
  probs[, "B"] <- apply(prob[,c("NDDB","NDEB"), drop = FALSE], 1, sum)
  probs[, "C"] <- apply(prob[,c("NDDC","NDEC","NDBC","NEBC"), drop = FALSE], 1, sum)

  if(detail)
    prob
  else
    probs
}
###############################################################################
calccounts <- function(rf = 0.5, gen = 10, trans = calctrans(rf, gen),
                       detail = FALSE)
{
  ## Expected counts.
  ## Nxyz = probability * count ending at z passing through x and later y.

  r2 = rf * rf;
  w2 = (1.0 - rf) * (1.0 - rf);
  rw = rf * (1.0 - rf);

  count <- matrix(0, gen, 12)
  dimnames(count) <- list(seq(gen),
                          c("NDDA","NDEA","NDBA","NEBA","NBA",
                            "NDDB","NDEB",
                            "NDDC","NDEC","NDBC","NEBC","NBC"))
  counts <- matrix(0, gen, 5)
  dimnames(counts) <- list(seq(gen), LETTERS[1:5])

  counts[, c("D","E")] <- trans[, c("PDD","PDE")] * trans[, c("MDD","MDE")]
  if(gen > 1) for(t in seq(2,gen)) {
    count[t,"NDDA"] <- sum(trans[seq(t-1),"PDD"] * (w2 / 4) *
                           trans[seq(t-1),"MDD"])
    count[t,"NDDB"] <- sum(trans[seq(t-1),"PDD"] * (rw / 2) * 2 ^ -seq(t-2,0) *
                           (1 + trans[seq(t-1),"MDD"]))
    count[t,"NDDC"] <- sum(trans[seq(t-1),"PDD"] * (r2 / 4) *
                           (2 + trans[seq(t-1),"MDD"]))
    if(t > 2) {
      count[t,"NDEA"] <- sum(trans[seq(t-1),"PDE"] * (r2 / 4) *
                             (2 + trans[seq(t-1),"MDE"]))
      count[t,"NDEB"] <- sum(trans[seq(t-1),"PDE"] * (rw / 2) * 2 ^ -seq(t-2,0) *
                             (1 + trans[seq(t-1),"MDE"]))
      count[t,"NDEC"] <- sum(trans[seq(t-1),"PDE"] * (w2 / 4) *
                             trans[seq(t-1),"MDE"])
      count[t,"NDBA"] <- 2 * sum(trans[seq(t-2),"PDD"] * (rw / 2) * trans[seq(t-2,1),"PBA"] *
                                 (1 + trans[seq(t-2),"MDD"]))
      count[t,"NDBC"] <- 2 * sum(trans[seq(t-2),"PDD"] * (rw / 2) * trans[seq(t-2,1),"PBA"] *
                                 (1 + trans[seq(t-2),"MDD"]))
      if(t > 3) {
        count[t,"NEBA"] <- 2 * sum(trans[seq(t-2),"PDE"] * (rw / 2) * trans[seq(t-2,1),"PBA"] *
                                   (1 + trans[seq(t-2),"MDE"]))
        count[t,"NEBC"] <- 2 * sum(trans[seq(t-2),"PDE"] * (rw / 2) * trans[seq(t-2,1),"PBA"] *
                                   (1 + trans[seq(t-2),"MDE"]))
      }
    }
  }

  count[, "NBA"] <- apply(count[,c("NDBA","NEBA"), drop = FALSE], 1, sum)
  count[, "NBC"] <- apply(count[,c("NDBC","NEBC"), drop = FALSE], 1, sum)

  counts[, "A"] <- apply(count[,c("NDDA","NDEA","NDBA","NEBA"), drop = FALSE], 1, sum)
  counts[, "B"] <- apply(count[,c("NDDB","NDEB"), drop = FALSE], 1, sum)
  counts[, "C"] <- apply(count[,c("NDDC","NDEC","NDBC","NEBC"), drop = FALSE], 1, sum)


  if(detail)
    count
  else
    counts
}

###############################################################################
###############################################################################
countit <- function(rf = 0.5, t = 3) ## works for t=3 and 4.
{
  raws <- countraw(rf, t)
  cook <- calccook(rf, t)

  list(rawp=raws$rawp,raw = raws$raw, cook = cook)
}
###############################################################################
countraw <- function(rf = 0.5, t = 3)
{
  rawp <- raw <- NA
  if(t == 3) {
    rawp <- c(dda = ((1-rf)^4 / 8 + ## DDA
                     (1-rf)^2 / 4), ## DAA
              dba = 2 * rf * (1 - rf) / 8, ## DBA (2 possible Bs for each A)
              dea = rf^4 / 8, ## DEA
              ddc = ((1-rf)^2 * rf^2 / 8 + ## DDC
                     rf^2 / 4), ## DCC
              dbc = 2 * rf * (1 - rf) / 8, ## DBC (2 possible Bs for each A)
              dec = rf^2 * (1-rf)^2 / 8) ## DEC
    raw <- rawp * c(0,1,4,2,1,2)
    rawp["tota"] <- sum(rawp[1:3])
    raw["tota"] <- sum(raw[1:3])
    rawp["totc"] <- sum(rawp[4:6])
    raw["totc"] <- sum(raw[4:6])
  }
  if(t == 4) {
    rawp <- c(dda = (((1-rf)^6 / 16) + ## DDDA
                     ((1-rf)^4 / 8) + ## DDAA
                     ((1-rf)^2 / 4) + ## DAAA
                     (rf^4 * (1-rf)^2 / 16)), ## DEDA
              dba = 2 * (((1-rf)^2 * rf * (1-rf) / 16) + ## DDBA
                (rf * (1-rf) / 16) + ## DBBA
                (rf * (1-rf) / 8) + ## DBAA
                (rf^2 * rf * (1-rf) / 16)), ## DEBA
              dea = ((rf^4 * (1-rf)^2 / 16) + ## DEEA
                     (rf^4 * (1-rf)^2 / 16) + ## DDEA
                     (rf^4  / 8)), ## DEAA
              ddc = (((1-rf)^4 * rf^2 / 16) + ## DDDC
                     ((1-rf)^2 * rf^2 / 8) + ## DDCC
                     (rf^2 / 4) + ## DCCC
                     (rf^6 / 16)), ## DEDC
              dbc = 2 * (((1-rf)^2 * rf * (1-rf) / 16) + ## DDBC
                (rf * (1-rf) / 16) + ## DBBC
                (rf * (1-rf) / 8) + ## DBCC
                (rf^2 * rf * (1-rf) / 16)), ## DEBC
              dec = ((rf^2 * (1-rf)^4 / 16) + ## DEEC
                     (rf^2 * (1-rf)^4 / 16) + ## DDEC
                     (rf^2 * (1-rf)^2 / 8))) ## DECC
    raw <- c(dda = 4 * rf^4 * (1-rf)^2 / 16, ## 0*(DDDA,DDAA,DAAA),4*DEDA
             dba = 2 * ((((1-rf)^2 * rf * (1-rf) / 16) + ## 1*DDBA
               (rf * (1-rf) / 16) + ## 1*DBBA
               (rf * (1-rf) / 8) + ## 1*DBAA
               3 * (rf^2 * rf * (1-rf) / 16))), ## 3*DEBA
             dea = 4 * rf^4 * ((2 * (1-rf)^2 / 16) + (1/8)), ## 4*(DEEA,DDEA,DEAA)
             ddc = (2 * (((1-rf)^4 * rf^2 / 16) + ## 2*DDDC
                         ((1-rf)^2 * rf^2 / 8) + ## 2*DDCC
                         (rf^2 / 4)) + ## 2*DCCC
                    6 * rf^6 / 16), ## 6*DEDA
             dbc = 2 * ((((1-rf)^2 * rf * (1-rf) / 16) + ## 1*DDBA
               (rf * (1-rf) / 16) + ## 1*DBBA
               (rf * (1-rf) / 8) + ## 1*DBAA
               3 * (rf^2 * rf * (1-rf) / 16))), ## 3*DEBA
             dec = 2 * ((rf^2 * (1-rf)^4 / 16) + ## 2*DEEC
               (rf^2 * (1-rf)^4 / 16) + ## 2*DDEC
               (rf^2 * (1-rf)^2 / 8))) ## 2*DECC
    rawp["tota"] <- sum(rawp[1:3])
    raw["tota"] <- sum(raw[1:3])
    rawp["totc"] <- sum(rawp[4:6])
    raw["totc"] <- sum(raw[4:6])
 }
 list(rawp=rawp, raw=raw)
}
###############################################################################
calccook <- function(rf, t)
{
  beta <- (rf^2 + (1-rf)^2)/2
  gamma <- ((1-rf)^2 - rf^2)/2
  alpha <- rf^2 / (2 * beta)
  SFt1 <- (1 - beta ^(t-2)) / (1 - beta)
  SF2t1 <- 2^(2-t) * (1 - (2*beta)^(t-2)) / (1-2*beta)

  k1b <- 0
  if(beta > 0)
    k1b <- kptothek(t-1, beta, beta^(t-1)) / beta
  k2b <- k2b2 <- 0
  if(beta > 0) {
    if(t > 2) {
      k2b <- kptothek(t-2, beta)
      k2b2 <- kptothek(t-2, 2 * beta)
    }
  }

  k1g <- 0
  SGt1 <- 0
  if(gamma > 0) {
    k1g <- kptothek(t-1, gamma, gamma^(t-1)) / gamma ^ 3
    SGt1 <- (1 - gamma^(t-2)) / (1-gamma)
  }

  ## D->..A
  dda <- 0
  if(t > 3)
    dda <- .125 * rf^2 * (1-rf)^2 * (k1b - 1 - k1g)

  dba1 <- rf * (1-rf) * (SFt1 - SF2t1) / 2
  dba2 <- 0.5 * rf * (1-rf) * (rf * (1-rf) + alpha) * k2b
  dba3 <- 0.5 * 2^(2-t) * rf * (1-rf) * alpha * k2b2
  dba <- dba1 + dba2 - dba3

  dea1 <- .125 * rf^4 * (k1b + 1 + k1g)
  dea2 <- 0.25 * rf^2 * (beta * SFt1 - gamma * SGt1)
  dea <- dea1 + dea2

  ## D->..C
  ddc1 <- 0
  #if(t > 3)
    ddc1 <- .125 * rf^2 * (rf)^2 * (k1b - 1 - k1g)
  ddc2 <- 0.25 * rf^2 * (beta * SFt1 + gamma * SGt1) ## should be prob(ddc) * 2 for t=4
  ddc <- ddc1 + ddc2

  dbc1 <- rf * (1-rf) * (SFt1 - SF2t1) / 2
  dbc2 <- 0.5 * rf * (1-rf) * (rf * (1-rf) + alpha) * k2b
  dbc3 <- 0.5 * 2^(2-t) * rf * (1-rf) * alpha * k2b2
  dbc <- dbc1 + dbc2 - dbc3

  dec1 <- 0
#  if(t > 3)
    dec1 <- .125 * rf^2 * (1-rf)^2 * (k1b + 1 + k1g)
  dec <- dec1

  c(dda=dda,ddc = ddc, dbc = dbc, dec = dec, totc = ddc+dbc+dec, ddc1=ddc1,ddc2=ddc2,
            k1b=k1b,k1g=k1g,SF2t1=SF2t1,SFt1=SFt1,SGt1=SGt1 #,dba1 = dba1, dba2 = dba2, dba3 = dba3
            )
}

