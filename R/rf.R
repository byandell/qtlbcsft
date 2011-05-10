######################################################################
#
# rf.R
#
# copyright (c) 2011, Brian S Yandell
# last modified May, 2011
# first written Apr, 2011
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
# 
# Part of the R/qtlbcsft package
# Contains: bcsft.rf,bcsft.rfrec,rfrec.plot,cmp.rf
#
######################################################################

######################################################################
## C computation of expected number of recombinations and log pr(gen1,gen2).
bcsft.rf <- function(rf = 0.5, s = 0, t = 2)
{
  n.gen <- 5
  if(t == 0)
    n.gen <- 2
  nrec2.ft <- nrec2.f2 <- matrix(NA, n.gen, n.gen)
  logprec.ft <- logprec.f2<- matrix(NA, n.gen, n.gen)
  for(gen1 in 1:n.gen) for(gen2 in 1:n.gen) {
    z <- myrf(gen1, gen2, rf, s, t)
    nrec2.ft[gen1,gen2] <- z[1]
    nrec2.f2[gen1,gen2] <- z[2]
    logprec.ft[gen1,gen2] <- z[3]
    logprec.f2[gen1,gen2] <- z[4]
  }
  list(nrec2.ft = nrec2.ft, nrec2.f2 = nrec2.f2, prec.ft = exp(logprec.ft), prec.f2 = exp(logprec.f2))
}

######################################################################
## Plot of rf estimated using expected count.
rfrec.plot <- function(s,t, add = FALSE, add.scheme = TRUE, ...)
{
  rf <- seq(0.01,0.49,by=0.02)
  rfrec <- apply(as.matrix(rf), 1, bcsft.rfrec,s,t)
  if(!add) {
    plot(c(0,.5),c(0,.5), type = "n", xlab = "recombination rate",
         ylab = "nrec2-based estimate")
    abline(0,1,col="gray")
  }
  lines(rf, rfrec, ...)
  if(add.scheme)
    text(0.49, max(rfrec), paste(s, t, sep = "."), adj = 0, cex = 0.65)
}

######################################################################
## Naive computation of rf using expected counts.

bcsft.rfrec <- function(rf = 0.5, s = 0, t = 2, n.ind = 1000,
                        ct = round(bcsft.prob(rf, s, t) * n.ind))
{
  ## Compute rf naively as in R/qtl's est.rf.
  ## This works exactly when only one recombination after F1 (e.g. BC1, F2).
  ## but not necessarily in general.
  
  n.gen <- 3
  if(t == 0)
    n.gen <- 2
  n.mei <- (t - (s == 0)) * 2 + s
  nrec2 <- matrix(NA, n.gen, n.gen)
  for(gen1 in 1:n.gen) for(gen2 in 1:n.gen) {
    z <- myrf(gen1, gen2, rf, s, t)
    nrec2[gen1,gen2] <- z[1]
  }
  sum(nrec2 * ct) / (n.mei * sum(ct))
}
bcsft.rfrecs <- function(s,t,ct, tol = 1e-6)
{
  rf <- bcsft.rfrec(0.5,s,t,ct = ct)
  repeat {
    new.rf <- bcsft.rfrec(rf,s,t,ct = ct)
    if(abs(rf - new.rf) < tol)
      break
    rf <- new.rf
  }
  rf
}
myrf <- function(gen1 = 1, gen2 = 1, rf = 0.5, s = 0, t = 2)
{
  cross.scheme <- c(s,t)
  ret <- double(4)
  z <- .C("rf_wrap", as.integer(gen1), as.integer(gen2), as.double(rf), as.integer(cross.scheme),
          ret = as.double(ret))
  z$ret
}

###########################################################################################
## Compare two rf objects. Messy, but useful.

cmp.rf <- function(rf1, rf2, tol = 1e-6)
{
  tmp <- row(rf1$rf) < col(rf1$rf)
  cat("Compare diagonals = number of meioses.\n")
  print(summary(c(diag(rf1$rf) - diag(rf2$rf))))
  cat("Compare upper triangles = LODs.\n")
  tmp2 <- rf1$rf[tmp] - rf2$rf[tmp]
  print(summary(tmp2)) ## upper tri: LOD
  tmp3 <- ((abs(tmp2) > tol) & !is.na(tmp2)) | ((is.na(rf1$rf[tmp]) + is.na(rf2$rf[tmp])) == 1)
  lod.diff <- NULL
  if(any(tmp3)) {
    rows <- row(tmp)[tmp][tmp3]
    cols <- col(tmp)[tmp][tmp3]
    lod.diff <- cbind(row=rows, col=cols, lod.1=rf1$rf[tmp][tmp3], lod.2=rf2$rf[tmp][tmp3],
                      rf.1=rf1$rf[t(tmp)][tmp3], rf.2=rf2$rf[t(tmp)][tmp3])
  }
  cat("Compare lower triangles = rfs.\n")
  tmp2 <- rf1$rf[t(tmp)] - rf2$rf[t(tmp)]
  print(summary(tmp2)) ## lower tri: rf
  tmp3 <- ((abs(tmp2) > tol) & !is.na(tmp2)) | ((is.na(rf1$rf[t(tmp)]) + is.na(rf2$rf[t(tmp)])) == 1)
  rf.diff <- NULL
  if(any(tmp3)) {
    rows <- row(tmp)[t(tmp)][tmp3]
    cols <- col(tmp)[t(tmp)][tmp3]
    rf.diff <- cbind(row=rows, col=cols, rf.1=rf1$rf[t(tmp)][tmp3], rf.2=rf2$rf[t(tmp)][tmp3],
                     lod.1=rf1$rf[tmp][tmp3], lod.2=rf2$rf[tmp][tmp3])
  }

  ## Are NAs for LOD and rf are different using flanking markers?
  cat("\nNAs for first: LOD vs. rf.")
  print(table(is.na(rf1$rf[tmp]), is.na(rf1$rf[t(tmp)])))
  cat("\nNAs for second: LOD vs. rf.")
  print(table(is.na(rf2$rf[tmp]), is.na(rf2$rf[t(tmp)])))

  miss.val <- function(rf1, rf2, tmp, type) {
    ## But some LODs that should be NA are actual values.
    tmp2 <- is.na(rf2$rf[tmp])
    if(any(tmp2)) {
      cat(paste(sum(tmp2), " ", type, "s that should be NA are actual values.\n", sep = ""))
      if(sum(tmp2) < 10)
        print(rf1$rf[tmp][tmp2])
      else
        print(summary(rf1$rf[tmp][tmp2]))
    }
  }
  miss.val(rf1, rf2, tmp, "LOD")
  miss.val(rf2, rf1, tmp, "LOD")
  miss.val(rf1, rf2, t(tmp), "rf")
  miss.val(rf2, rf1, t(tmp), "rf")

  cat("\nNA LODs vs. 0 rfs. for first.")
  print(table(is.na(rf1$rf[tmp]), is.na(rf1$rf[t(tmp)]) | (rf1$rf[t(tmp)] %in% c(0,1))))
  cat("\nNA LODs vs. 0 rfs. for second.")
  print(table(is.na(rf2$rf[tmp]), is.na(rf2$rf[t(tmp)]) | (rf2$rf[t(tmp)] %in% c(0,1))))

  cat("\nCompare NA for LODs: first vs. second.")
  print(table(is.na(rf1$rf[tmp]), is.na(rf2$rf[tmp])))
  cat("\nCompare NA vs. 0 for LODs: first vs. second.")
  print(table(is.na(rf1$rf[tmp]), rf2$rf[tmp] == 0))
  cat("\nCompare NA or 0 vs. 0 for LODs: first vs. second.")
  print(table(is.na(rf1$rf[tmp]) + (!is.na(rf1$rf[tmp]) & (rf1$rf[tmp] == 0)), rf2$rf[tmp] == 0))

  list(lod.diff = lod.diff, rf.diff = rf.diff)
}
