## r 0.5 0.5
## t 2   3
## A 0   2/3
## B 1   4/3
## C 2   14/9
## D 0   2
## E 2   2

## For compatibility with C.
R_pow <- function(a, b) a ^ b
#################################################################################################
## Not used.
expect_bcs <- function(rf = 0.5, s = 2, c = count_bcs(rf, s), p = prob_bcs(rf, s))
{
  out <- c / p
  out[is.na(out)] <- 0
  out
}
expect_ft <- function(rf = 0.5, t = 2, c = count_ft(rf, t), p = prob_ft(rf, t))
{
  out <- c / p
  out[is.na(out)] <- 0
  out
}
probs <- function(rf = 0.5, bc.gen = 2, gen = 10)
{
  out <- matrix(0, bc.gen + gen, 7)
  if(bc.gen > 0) {
    for(s in seq(bc.gen))
      out[s,] <- tmp <- prob_bcsft(rf, s, 0)
    for(t in seq(gen))
      out[bc.gen + t,] <- tmp <- prob_bcsft(rf, bc.gen, t)
  }
  else {
    for(t in seq(gen))
      out[t,] <- tmp <- prob_bcsft(rf, 0, t)
  }
  dimnames(out) <- list(gennames(c(bc.gen,gen)), names(tmp))
  out
}
counts <- function(rf = 0.5, bc.gen = 2, gen = 10)
{
  out <- matrix(0, bc.gen + gen, 7)
  if(bc.gen > 0) {
    for(s in seq(bc.gen))
      out[s,] <- tmp <- count_bcsft(rf, s, 0)
    for(t in seq(gen))
      out[bc.gen + t,] <- tmp <- count_bcsft(rf, bc.gen, t)
  }
  else {
    for(t in seq(gen))
      out[t,] <- tmp <- count_bcsft(rf, 0, t)
  }
  dimnames(out) <- list(gennames(c(bc.gen,gen)), names(tmp))
  out
}
expects <- function(rf = 0.5, bc.gen = 2, gen = 10)
{
  out <- matrix(0, bc.gen + gen, 7)
  if(bc.gen > 0) {
    for(s in seq(bc.gen))
      out[s,] <- tmp <- expect_bcsft(rf, s, 0)
    for(t in seq(gen))
      out[bc.gen + t,] <- tmp <- expect_bcsft(rf, bc.gen, t)
  }
  else {
    for(t in seq(gen))
      out[t,] <- tmp <- expect_bcsft(rf, 0, t)
  }
  dimnames(out) <- list(gennames(c(bc.gen,gen)), names(tmp))
  out
}

############################################################################################
############################################################################################
prob_bcs <- function(rf = 0.5, s = 2)
{
  #######BCs probabilities:
  if(s > 0) {
    w = 1 - rf
    ws = R_pow(w,s);
    s2 = R_pow(2,s);
    s1 = s - 1;

    ## Ls + 2*Ms + Ns = 1.0 regardless of s
    PbD = ws / s2;
    PbB = (1 - ws)/ s2;
    PbA = (s2 - 2 + ws)/ s2;
    
    c(A1=PbA, B1=PbB, C=0, D=PbD, E=0, A0=0, B0=0)
  }
  else
    c(A1=0, B1=0, C=0, D=1, E=0, A0=0, B0=0)
}
prob_ft <- function(rf = 0.5, t = 2)
{
  ## compute transition probabilities to leave double het states
  t1 = t - 1.0;
  t2 = R_pow(2.0, t); ## 2^t
  w = 1.0 - rf;
  w2 = w * w;
  r2 = rf * rf;
  rw = w * rf;
  w2pr2 = (w2 + r2);

  ## transient state Bt
  Bt = (2.0 * rw / t2) * ((1.0 - R_pow(w2pr2, t1)) / (1 - w2pr2));

  ## A = prob go from AB.ab to AB.AB at step t
  ## D1 -> Dk or Ek -> Ak+1 -> At OR D1 -> Dk or Ek -> Bk+1 -> Ak+s -> At
  beta = w2pr2 / 2.0; 
  gamma = (w2 - r2) / 2.0; 
  beta1 = R_pow(beta, t1);
  ## Ft = prob still in D or E at step t
  ## SFt = sum of Fk from 1 to t-1
  ## Dt and Et depend on Ft and Gt
  SFt = (1.0 - beta1) / (1.0 - beta);
  SGt = (1.0 - R_pow(gamma, t1)) / (1.0 - gamma);
  SDt = (SFt + SGt) / 2.0;
  SEt = (SFt - SGt) / 2.0;
  ## sum Fk * rw/2 * prob(B->A in remaining steps)
  SFta = (rw / 2.0) * (SFt - 2.0 * ((2.0 / t2) - beta1) / (1.0 - 2.0 * beta));
  ## note symmetry in At and Ct
  At = SDt * w2 / 4 + SEt * r2 / 4 + SFta;
  Ct = SDt * r2 / 4 + SEt * w2 / 4 + SFta;

  ## compute transition probabilities to stay in double het states
  t1 = t - 1.0;
  w = 1.0 - rf;
  w2 = w * w;
  r2 = rf * rf;

  Ft = R_pow((w2 + r2) / 2.0, t1);
  Gt = R_pow((w2 - r2) / 2.0, t1);

  ## now the terms we really want!
  Dt = (Ft + Gt) / 2.0;
  Et = (Ft - Gt) / 2.0;

  ## Both A's have same probability for Ft.
  c(A1=At, B1=Bt, C=Ct, D=Dt, E=Et, A0=At, B0=Bt)
}
prob_bcsft <- function(rf = 0.5, s = 2, t = 2,
                       p.bcs = prob_bcs(rf, s), 
                       p.ft = prob_ft(rf, t + (s > 0)))
{
  if(s == 0)
    return(p.ft)
  if(t == 0)
    return(p.bcs)

  ## BCsFt (watch out for t=1 situation *****)
  ## probabilities: (Pxx)
  ## D->x: PfDx using calculations from before, but need to multiple by PbD = BCs probability
  ## 4 states now split into 7 (see above)

  ## PbDfx = PbD * PfDx (PfDx = probability of D->x in Ft)

  PbDfD = p.bcs["D"] * p.ft["D"];
  PbDfE = p.bcs["D"] * p.ft["E"];
  PfD = PbDfD;
  PfE = PbDfE;

  ## B: distinguish B1 (reachable in BCs and Ft) and B0 (reachable only in Ft)
  ## (watch out for t=1 situation)
  ## add in D->B calcs for Ft

  PbDfB = p.bcs["D"] * p.ft["B0"]
  PBfB = R_pow(0.5, t)
  PbBfB1 = p.bcs["B1"] * PBfB; 
  PfB0 = PbDfB;
  PfB1 = PbDfB + PbBfB1;

  ## B->C: go to B in BCs, then later migrate to C
  ## (there are 2 separate C's to consider, one for each B)

  PbDfC = p.bcs["D"] * p.ft["C"]
  ## PBfC = sum from k=0 to t-1 (1/2)^k * (1/4)
  PBfC = 0.25 * sum(R_pow(0.5, seq(0, t-1)));
  PbBfC = p.bcs["B1"] * PBfC;
  PfC = PbDfC + PbBfC;

  ## B->A: same as B->C but * 2 (both B's go to the same A:ab.ab)
  PbBfA1 = 2 * PbBfC;

  ## A: ab.ab
  PbAfA1 = p.bcs["A1"] * 1; ## (no further change)

  PbDfA = p.bcs["D"] * p.ft["A0"]
  PfA0 = PbDfA;
  PfA1 = PbDfA + PbBfA1 + PbAfA1;

  out <- c(A1=PfA1, B1=PfB1, C=PfC, D=PfD, E=PfE, A0=PfA0, B0=PfB0)
  names(out) <- c("A1","B1","C","D","E","A0","B0")
  out
}
###########################################################################################
count_bcs <- function(rf = 0.5, s = 2, prob = prob_bcs(rf, s))
{
  transgen = rep(0, 7);
  names(transgen) <- c("A1","B1","C","D","E","A0","B0")

  PbB = prob["B1"]
  PbD = prob["D"]

  ## PbA = PbA1 + PbBA1;
  PbA1 = (1 - rf) * (1 - PbD) / (1 + rf);
  ##  2 ways to get to A1 from B:Ab.ab or B:aB.ab
  PbBA1 = 1 - PbD - 2 * PbB - PbA1 ;    ## stay in D for 0 to k generations, then B for 1 to j generations, then A

  transgen[1+1] = PbB; ## D->B
  transgen[1+0] = PbBA1; ## D->A
  
  transgen;
}

count_ft <- function(rf = 0.5, t = 2)
{
  transgen = rep(0, 7);
  names(transgen) <- c("A1","B1","C","D","E","A0","B0")
  if(t < 2)
    return(transgen)

  t1 = t - 1.0;
  t2 = R_pow(2.0, -t1);

  r2 = rf * rf;
  w2 = (1.0 - rf) * (1.0 - rf);
  rw = rf * (1.0 - rf);

  alpha = r2 / (w2 + r2);

  ## beta = probability of D or E at each step.
  beta = (w2 + r2) / 2.0;
  beta1 = R_pow(beta, t1);
  beta2 = 1.0;
  if(t > 2.0) beta2 = beta1 / beta;
  ## Calculations for D->DE->B->AC
  sbeta1 = (1.0 - beta1) / (1.0 - beta); ## SFt
  sbeta2 = 0.0;
  if(t > 2.0) sbeta2 = (1.0 - beta1 / beta) / (1.0 - beta); ## SF(t-1)
  s2beta1 = (t2 - beta1) / (1.0 - 2.0 * beta); ## sum from 1 to t-1 of of (2*beta)^(k-1).
  s2beta2 = (2.0 * t2 - (beta1 / beta)) / (1.0 - 2.0 * beta); ## sum from 1 to t-2 of of (2*beta)^(k-1).
  
  
  gamma = (w2 - r2) / 2.0; 
  gamma1 = 1.0;
  if(t > 1) gamma1 = R_pow(gamma, t1);
  gamma2 = 1.0;
  if(t > 2) gamma2 = R_pow(gamma, t1 - 1);
  sgamma1 = 1;
  sgamma2 = 0;
  s2gamma1 = t2;
  s2gamma2 = 0;
  if(t > 1) {
    sgamma2 = 1;
    s2gamma2 = t2 * 2.0;
  }
  if(gamma > 0) {
    sgamma1 = (1.0 - gamma1) / (1.0 - gamma); ## SGt
    sgamma2 = (1.0 - gamma2) / (1.0 - gamma); ## SG(t-1)
    s2gamma1 = (t2 - gamma1) / (1.0 - 2.0 * gamma); ## sum from 1 to t-1 of of (2*gamma)^(k-1).
    s2gamma2 = (2.0 * t2 - (gamma1 / gamma)) / (1.0 - 2.0 * gamma); ## sum from 1 to t-2 of of (2*gamma)^(k-1).
  }

  ## kptothek(t,p) = sum of k * p^k from 1 to t-1.
  k1b = kptothek(t1, beta, beta1) / beta;
  k2b = t2 * kptothek(t1, 2.0 * beta, beta1 / t2) / (2 * beta);
  k1g = 0.0;
  k2g = 0.0;
  k1g2 = 0.0;
  k2g2 = 0.0;
  k1b2 = 0.0;
  k2b2 = 0.0;
  if(t > 2) {
    k1g = 1.0;
    k2g = t2;
    if(t > 3) {
      k1g2 = 1.0;
      k2g2 = 2.0 * t2;
    }
    k1b2 = kptothek(t1 - 1.0, beta, beta2) / beta;
    k2b2 = 2.0 * t2 * kptothek(t1 - 1.0, 2.0 * beta, beta2 / (2 * t2)) / (2 * beta);
  }
  if(gamma > 0) {
    ## Possible savings in doing the sum...
    k1g = kptothek(t1, gamma, gamma1) / gamma;
    k2g = t2 * kptothek(t1, 2.0 * gamma, gamma1 / t2) / (2.0 * gamma);
    k1g2 = kptothek(t1 - 1.0, gamma, gamma2) / gamma;
    k2g2 = 2.0 * t2 * kptothek(t1 - 1.0, 2.0 * gamma, gamma2 / (2.0 * t2)) / (2.0 * gamma);
  }

  ## 0: AB.AB or ab.ab (At)
  ## 1: Ab.Ab or aB.aB (Ct)
  ## 2: rest (8 of them) (Bt)
  ## 3: AB.ab or ab.AB (Dt)
  ## 4: Ab.aB or aB.Ab (Et)

  ## Below, x->y->z->w refers to generations 1, k, k+1, t
  ## In case of state at k+1 is B, there is a further transition to A or C at step k+s summed out

  ## D: based on odd counts
  ## E: based on even counts
  transgen[1+1] = rw * (s2beta1 +  2.0 * r2 * k2b); ## D->DE->B->B1
  transgen[1+3] = 0.5 * t1 * r2 * (beta2 - gamma2); ## D->D
  transgen[1+4] = 0.5 * t1 * r2 * (beta2 + gamma2); ## D->E

  ## End point is A or C. Newly revised.
  ndda = (r2 / 2) * (k1b - k1g);
  NDDA = (w2 / 4) * ndda;
  NDDC = (r2 / 4) * (ndda + (sbeta1 + sgamma1));

  NDEA = 0.0;
  NDEC = 0.0;
  NDBA = 0.0;
  NEBA = 0.0;
  if(t > 2.0) {
    ndea = (r2 / 2.0) * (k1b + k1g);
    NDEA = (r2 / 4.0) * (ndea + (sbeta1 - sgamma1));
    NDEC = (w2 / 4.0) * ndea;
    
    nbab = rw * (0.25 * (sbeta2 - s2beta2) + (r2 / 2.0) * (0.5 * k1b2 - k2b2));
    nbag = rw * (0.25 * (sgamma2 - s2gamma2) - (r2 / 2.0) * (0.5 * k1g2 - k2g2));
    NDBA = nbab + nbag;
    if(t > 3.0)
      NEBA = nbab - nbag;
  }
  transgen[1+0] = (NDDA + NDEA + NDBA + NEBA); ## D->A1
  transgen[1+2] = (NDDC + NDEC + NDBA + NEBA); ## D->C

  transgen[1+5] = transgen[1+0]; ## D->A0
  transgen[1+6] = transgen[1+1]; ## D->B0
  transgen;
}
count_bcsft <- function(rf = 0.5, s = 2, t = 2,
                        p.bcs = prob_bcs(rf, s), 
                        n.bcs = count_bcs(rf, s), 
                        n.ft = count_ft(rf, t+(s>0)))
{
  transgen = rep(0, 7);
  names(transgen) <- c("A1","B1","C","D","E","A0","B0")

  if(s == 0)
    return(n.ft)
  if(t == 0)
    return(n.bcs)

  ## numerator of counts:
  ## NbDfx = PbD * NfDx (NfDx = numerator of cound for D->x in Ft)

  ## B: 
  ## NbDfB = NbDfB0 = NbDfB1 = PbD * NfDB (counts as before)
  NfB0 = p.bcs["D"] * n.ft["B0"];
  NfB1 = NfB0 + p.bcs["B1"] * R_pow(0.5, t);

  ## C: 
  NbBfC = p.bcs["B1"] * 0.5 * (1 - R_pow(2, -t));
  NfC = p.bcs["D"] * n.ft["C"] + NbBfC;

  ## A:
  NfA0 = p.bcs["D"] * n.ft["A0"];
  NfA1 = NfA0 + 2 * NbBfC + n.bcs["A1"];

  ## Only As and C are correct!
  transgen[1+0] = NfA1; ## D->A1
  transgen[1+1] = NfB1; ## D->B1
  transgen[1+2] = NfC; ## D->C
  transgen[1+3] = p.bcs["D"] * n.ft["D"]; ## D->D
  transgen[1+4] = p.bcs["D"] * n.ft["E"]; ## D->E
  transgen[1+5] = NfA0; ## D->A0
  transgen[1+6] = NfB0; ## D->B0
  
  transgen;
}
expect_bcsft <- function(rf = 0.5, s = 2, t = 10, c = count_bcsft(rf, s, t),
                         p = prob_bcsft(rf, s, t))
{
  out <- c / p
  out[is.na(out)] <- 0
  out
}
#################################################################################################
kptothek <- function(t, p, ptothet = p ^ t)
{
  ## sum of k * p^k from 0 to t-1
  tmp = (1.0 - p);
  ((p - t * ptothet + (t - 1.0) * p * ptothet) / (tmp * tmp));
}
