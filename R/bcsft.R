##############################################################################
### Reachable in BCs:
### D : AB.ab
### B1: Ab.ab, aB.ab
### A1: ab.ab
### Reachable only in Ft:
### A0: AB.AB
### B0: AB.Ab, AB.aB
### C : Ab.Ab, aB.aB
### E : Ab.aB
##############################################################################
R_pow <- function(a, b) a ^ b

prob.bcs <- function(rf = 0.5, s = 2)
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
    ## PbA = PbA1 + PbBA1;
    PbA1 = w * (1 - PbD) / (2 - w);
    ##  2 ways to get to A1 from B:Ab.ab or B:aB.ab
    PbBA1 = 1 - PbD - 2 * PbB - PbA1 ;    ## stay in D for 0 to k generations, then B for 1 to j generations, then A
    c(A=PbA, B=PbB, D=PbD, A1 = PbA1, BA1 = PbBA1)
  }
  else
    c(A=0, B=0, D=1, A1=0, BA1=0)
}
probs.bcs <- function(rf = 0.5, bc.gen = 10)
{
  if(bc.gen < 1)
    stop("bc.gen must be at least 1")

  out <- matrix(0, bc.gen, 5)
  for(s in seq(bc.gen))
    out[s,] <- tmp <- prob.bcs(rf, s)
  dimnames(out) <- list(seq(bc.gen), names(tmp))
  out
}
#######################################################################
  
count.bcs <- function(rf = 0.5, s = 2, prob = prob.bcs(rf, s))
{
  if(s > 0) {
    PbA = prob["A"]
    PbB = prob["B"]
    PbD = prob["D"]
    PbA1 = prob["A1"]
    PbBA1 = prob["BA1"]
    ## PbA = PbA1 + PbBA1;
    PbA1 = (1 - rf) * (1 - PbD) / (1 + rf);
    ##  2 ways to get to A1 from B:Ab.ab or B:aB.ab
    PbBA1 = 1 - PbD - 2 * PbB - PbA1 ;    ## stay in D for 0 to k generations, then B for 1 to j generations, then A

    ## counts:
    count_D = 0;
    NbD = 0 * PbD;
    count_B = 1;
    NbB = 1 * PbB;

    ## have to count these separately:
    NbA1 = 0 * PbA1 
    NbBA1 = 1 * PbBA1;
    NbA = NbA1 + NbBA1; ## = 2 * PbBA1

    out <- c(NbA, NbB, NbD, NbA1, NbBA1)
  }
  else
    out <- rep(0, 5)
  names(out) <- c("A","B","D","A1","BA1")
  out
}
counts.bcs <- function(rf = 0.5, bc.gen)
{
  if(bc.gen < 1)
    stop("bc.gen must be at least 1")

  out <- matrix(0, bc.gen, 5)
  for(s in seq(bc.gen))
    out[s,] <- tmp <- count.bcs(rf, s)
  dimnames(out) <- list(seq(bc.gen), names(tmp))
  out
}
expect.bcs <- function(rf = 0.5, s = 2, prob = prob.bcs(rf, s),
                       count = count.bcs(rf, s))
{
  ## expected counts = NbA/PbA = 2 * PbBA1 / (PbA)
  ## = 2 * (1 - PbD - 2 * PbB - PbA1) / (1 - PbD - 2 * PbB)
  ## = 2 - PbA1 / (1 - PbD - 2 * PbB)

  out <- c(count / prob, total = sum(count[1:3] / prob[1:3]))
  out[is.na(out)] <- 0
  out
}
expects.bcs <- function(rf = 0.5, bc.gen)
{
  out <- matrix(0, bc.gen, 6)
  for(s in seq(bc.gen))
    out[s,] <- tmp <- expect.bcs(rf, s)
  dimnames(out) <- list(seq(bc.gen), names(tmp))
  out
}
##############################################################################
prob.bcsft <- function(rf = 0.5, s = 2, t = 2,
                       p.bcs = prob.bcs(rf, s), 
                       p.ft = calcprobs(rf, t+(s>0)))
{
  extras <- c("A0", "B0")
  if(s == 0) {
    out <- p.ft[t,]
    names(out)[1:2] <- paste(LETTERS[1:2], 1, sep = "")
    seven <- rep(0, 2)
    names(seven) <- extras
    return(c(out, seven))
  }
  if(t == 0) {
    out <- rep(0, 7)
    names(out) <- c(LETTERS[1:5], extras)
    names(out)[1:2] <- paste(LETTERS[1:2], 1, sep = "")
    out[c(1,2,4)] <- p.bcs[c(1,2,3)]
    return(out)
  }

  ## BCsFt (watch out for t=1 situation *****)
  ## probabilities: (Pxx)
  ## D->x: PfDx using calculations from before, but need to multiple by PbD = BCs probability
  ## 4 states now split into 7 (see above)

  ## PbDfx = PbD * PfDx (PfDx = probability of D->x in Ft)

  PbDfD = p.bcs["D"] * p.ft[t+1, "D"];
  PbDfE = p.bcs["D"] * p.ft[t+1, "E"];
  PfD = PbDfD;
  PfE = PbDfE;

  ## B: distinguish B1 (reachable in BCs and Ft) and B0 (reachable only in Ft)
  ## (watch out for t=1 situation)
  ## add in D->B calcs for Ft

  PbDfB = p.bcs["D"] * p.ft[t+1, "B"]
  PBfB = R_pow(0.5, t)
  PbBfB1 = p.bcs["B"] * PBfB; 
  PfB0 = PbDfB;
  PfB1 = PbDfB + PbBfB1;

  ## B->C: go to B in BCs, then later migrate to C
  ## (there are 2 separate C's to consider, one for each B)

  PbDfC = p.bcs["D"] * p.ft[t+1, "C"]
  ## PBfC = sum from k=0 to t-1 (1/2)^k * (1/4)
  PBfC = 0.25 * sum(R_pow(0.5, seq(0, t-1)));
  PbBfC = p.bcs["B"] * PBfC;
  PfC = PbDfC + PbBfC;

  ## B->A: same as B->C but * 2 (both B's go to the same A:ab.ab)
  PbBfA1 = 2 * PbBfC;

  ## A: ab.ab
  PbAfA1 = p.bcs["A"] * 1; ## (no further change)

  PbDfA = p.bcs["D"] * p.ft[t+1, "A"]
  PfA0 = PbDfA;
  PfA1 = PbDfA + PbBfA1 + PbAfA1;

  five = c(PfA1, PfB1, PfC, PfD, PfE)
  names(five) <- paste(LETTERS[1:5], c(rep(1,2),rep("",3)), sep = "")
  seven <- c(PfA0, PfB0)
  names(seven) <- c("A0", "B0")
  c(five, seven) 
}
probs.bcsft <- function(rf = 0.5, bc.gen = 2, gen = 10)
{
  out <- matrix(0, bc.gen + gen, 7)
  if(bc.gen > 0) {
    for(s in seq(bc.gen))
      out[s,] <- tmp <- prob.bcsft(rf, s, 0)
    for(t in seq(gen))
      out[bc.gen + t,] <- tmp <- prob.bcsft(rf, bc.gen, t)
  }
  else {
    for(t in seq(gen))
      out[t,] <- tmp <- prob.bcsft(rf, 0, t)
  }
  dimnames(out) <- list(gennames(c(bc.gen,gen)), names(tmp))
  out
}

count.bcsft <- function(rf = 0.5, s = 2, t = 2,
                       p.bcs = prob.bcs(rf, s), 
                       n.bcs = count.bcs(rf, s), 
                       p.ft = calcprobs(rf, t+(s>0)),
                       n.ft = calccounts(rf, t+(s>0)))
{
  extras <- c("A0", "B0")
  if(s == 0) {
    out <- n.ft[t,]
    names(out)[1:2] <- paste(LETTERS[1:2], 1, sep = "")
    seven <- rep(0, 2)
    names(seven) <- extras
    return(c(out, seven))
  }
  if(t == 0) {
    out <- rep(0, 7)
    names(out) <- c(LETTERS[1:5], extras)
    names(out)[1:2] <- paste(LETTERS[1:2], 1, sep = "")
    out[c(1,2,4)] <- n.bcs[c(1,2,3)]
    return(out)
  }

  ## numerator of counts:
  ## D: 0 + Ls * (counts from before)
  ## NbDfx = PbD * NfDx (NfDx = numerator of cound for D->x in Ft)

  NbDfD = p.bcs["D"] * n.ft[t+1, "D"];
  NbDfE = p.bcs["D"] * n.ft[t+1, "E"];
  NfD = NbDfD;
  NfE = NbDfE;

  ## B: 
  ## NbDfB = NbDfB0 = NbDfB1 = PbD * NfDB (counts as before)
  NbDfB = p.bcs["D"] * n.ft[t+1, "B"];
  PBfB = R_pow(0.5, t)
  PbBfB1 = p.bcs["B"] * PBfB; 
  NbBfB1 = PbBfB1 * 1;
  NfB0 = NbDfB;
  NfB1 = NbDfB + NbBfB1;
  NfB = NfB0 + NfB1;

  ## C: 
  NbDfC = p.bcs["D"] * n.ft[t+1, "C"];
  PBfC = 0.25 * sum(R_pow(0.5, seq(0, t-1)));
  PbBfC = p.bcs["B"] * PBfC;
  NbBfC = PbBfC; 
  NfC = NbDfC + NbBfC;

  ## A:
  NbDfA = p.bcs["D"] * n.ft[t+1, "A"];
  NfA0 = NbDfA;
  NbAfA1 = n.bcs["A"];
  NbBfA1 = 2 * NbBfC;
  NfA1 = NbDfA + NbBfA1 + NbAfA1
  NfA = NfA0 + NfA1;

  five = c(NfA1, NfB1, NfC, NfD, NfE)
  names(five) <- paste(LETTERS[1:5], c(rep(1,2),rep("",3)), sep = "")
  seven <- c(NfA0, NfB0)
  names(seven) <- c("A0", "B0")
  c(five, seven) 
}
counts.bcsft <- function(rf = 0.5, bc.gen = 2, gen = 10)
{
  out <- matrix(0, bc.gen + gen, 7)
  if(bc.gen > 0) {
    for(s in seq(bc.gen))
      out[s,] <- tmp <- count.bcsft(rf, s, 0)
    for(t in seq(gen))
      out[bc.gen + t,] <- tmp <- count.bcsft(rf, bc.gen, t)
  }
  else {
    for(t in seq(gen))
      out[t,] <- tmp <- count.bcsft(rf, 0, t)
  }
  dimnames(out) <- list(gennames(c(bc.gen,gen)), names(tmp))
  out
}
#################
## Main work to do: 
## 1. find out NAA1
## 2. expand from 4 to 7 states
## 3. build new code for above
## 4. merge with existing code

checkprob <- function(rf = 0.5, s = 2, t = 10,
                      matp = genprobs(rf, cross.scheme = c(s,t)),
                      bcsft = probs.bcsft(rf, bc.gen = s, gen = t),
                      pr = probs(rf, s, t), roundoff = TRUE)
{
  tmp <- bcsft
  tmp[,"A0"] <- matp[,"A:AB.AB"]
  tmp[,"A1"] <- matp[,"A:ab.ab"]
  tmp[,"B0"] <- matp[,"B:AB.Ab"]
  tmp[,"B1"] <- matp[,"B:Ab.ab"]
  tmp[,c("C","D","E")] <- matp[,c("C:Ab.Ab","D:AB.ab","E:Ab.aB")]
  matp <- tmp

  tmp <- bcsft
  ## Get As.
  tmp[,"A0"] <- pr[,"A0"]
  tmp[,"A1"] <- pr[,"A1"]
  tmp[,"B0"] <- pr[,"B0"]
  tmp[,"B1"] <- pr[,"B1"]
  tmp[,c("C","D","E")] <- pr[,c("C","D","E")]
  
  out <- list(bcsft = bcsft, matrix = matp, qtlt = tmp, diff = matp - tmp)
  if(roundoff)
    out <- lapply(out, round, 6)
  out
}
checkcount <- function(rf = 0.5, s = 2, t = 10,
                       matp = sumcts(rf, cross.scheme = c(s,t), numerator = TRUE),
                       bcsft = counts.bcsft(rf, bc.gen = s, gen = t),
                       pr = counts(rf, s, t), roundoff=TRUE)
{
  checkprob(rf,s,t,matp,bcsft,pr, roundoff)
}
checkexp <- function(rf = 0.5, s = 2, t = 10, roundoff=TRUE)
{
  prob <- checkprob(rf,s,t, roundoff = FALSE)
  count <- checkcount(rf,s,t, roundoff = FALSE)
  for(i in c("bcsft","matrix")) {
    count[[i]] <- count[[i]] / prob[[i]]
    count[[i]][is.na(count[[i]])] <- 0
  }
  count$diff <- NULL
  tmp <- counts(rf, s, t)
  tmp <- prob$bcsft
  pr <- expects(rf, s, t)
  tmp[,"A0"] <- pr[,"A0"]
  tmp[,"A1"] <- pr[,"A1"]
  tmp[,"B0"] <- pr[,"B0"]
  tmp[,"B1"] <- pr[,"B1"]
  tmp[,c("C","D","E")] <- pr[,c("C","D","E")]
  count$qtlt <- tmp
  
  count$diff <- count$qtlt - count$matrix
  if(roundoff)
    count <- lapply(count, round, 6)
  count
}
