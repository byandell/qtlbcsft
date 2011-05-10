pr3 <- function(rf)
{
  p22 <- ((1 - rf) / 2) ^ 2
  p12 <- rf * (2 - rf) / 4
  matrix(c(0.5 + p22, rep(p12, 2), p22), 2, 2)
}

like3 <- function(rf, ct, s, t, phase = FALSE)
{
  pr <- bcsft.prob(rf, s, t, phase)
  sum(ct * log(pr))
}
bcsft.like <- function(rf = 0.2, s = 2, t = 0, n.ind = 1000,
                       ct = round(bcsft.prob(rf, s, t, phase) * n.ind),
                       phase = FALSE, add = FALSE, ylim = c(-0.1,0), ...)
{
  rr <- seq(0.01, 0.99, by = 0.01)
  ll <- apply(as.matrix(rr), 1, like3, ct, s, t, phase)
  if(!add) {
    plot(c(0,.5), c(-1,0), type = "n", xlab = "recombination rate", ylab = "(loglik-max)/range", ylim = ylim)
  }
  lines(rr, (ll - max(ll)) / diff(range(ll)), ...)
  wh <- which.max(ll)
  abline(v = rr[wh])
  c(rf, rr[wh], ll[wh], range(ll))
}
golden <- function(ct, fn = like3, cross.scheme, tol = 1e-6, ...)
{
  ## en.wikipedia.org/wiki/Golden_section_search

  s <- cross.scheme[1]
  t <- cross.scheme[2]

  phi <- (1 + sqrt(5)) / 2
  resphi <- 2 - phi

  x <- c(0, 1)
  y <- apply(as.matrix(x), 1, fn, ct, s, t, ...)
  wh <- which.max(y)
  newx <- x[wh] + resphi * (x[3 - wh] - x[wh])

  x <- c(x[wh], newx, x[3 - wh])
  y <- c(y[wh], fn(newx, ct, s, t, ...), y[3 - wh])

  outx <- x[c(1,3,2)]
  outy <- y[c(1,3,2)]

  ## x1 and x3 are the current bounds; the minimum is between them.
  ## x2 is the center point, which is closer to x1 than to x3

  repeat {
    ## Create a new possible center in the area between x2 and x3, closer to x2
    x[4] <- x[2] + resphi * (x[3] - x[2])

    ## Evaluate termination criterion
    if(abs(x[3] - x[1]) < tol)
      break
 
    y[4] <- fn(x[4], ct, s, t, ...)
    outx <- c(outx, x[4])
    outy <- c(outy, y[4])

    if(y[4] > y[2]) {
      x <- x[c(2,4,3)]
      y <- y[c(2,4,3)]
    }
    else {
      x <- x[c(4,2,1)]
      y <- y[c(4,2,1)]
    }
  }
  list(max = (x[3] + x[1]) / 2, x = outx, y = outy, ct = ct, cross.scheme = cross.scheme)
}

like.plot <- function(object, ...)
{
  cross.scheme = object$cross.scheme
  
  out <- bcsft.like(rf = 0.5, s = cross.scheme[1], t = cross.scheme[2], 
                    ct = object$ct, ...)
  range.y <- out[4:5]
  n.out <- length(object$x)
  abline(v = object$x, col = "gray")
  for(i in seq(n.out))
    points(object$x[i], (object$y[i] - range.y[2]) / diff(range.y))
  abline(v = object$max, col = "red")
  invisible(out)
}
