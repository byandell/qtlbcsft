M_LN2 <- log(2)
prob.ft <- function(rf, t, transpr)
{
  ## compute transition probabilities to leave double het states */
  t1 = t - 1.0;
  t2 = R_pow(2.0, t); ## 2^t */
  t1m2 = 2.0 / t2;
  w = 1.0 - rf;
  w2 = w * w;
  r2 = rf * rf;
  rw = w * rf;

  transpr <- rep(0.0, 10);

  ## A = prob go from AB.ab to AB.AB at step t */
  ## D1 -> Dk or Ek -> Ak+1 -> At OR D1 -> Dk or Ek -> Bk+1 -> Ak+s -> At */
  beta = (w2 + r2) / 2.0; 
  gamma = (w2 - r2) / 2.0; 
  beta1 = R_pow(beta, t1);
  gamma1 = R_pow((w2 - r2) / 2.0, t1);
  ## beta1 = prob still in D or E at step t */
  ## sbeta1 = sum of beta^(k-1) from 1 to t-1 */
  ## Dt and Et depend on sbeta1 and sgamma1 */
  sbeta1 = (1.0 - beta1) / (1.0 - beta);
  sgamma1 = (1.0 - R_pow(gamma, t1)) / (1.0 - gamma);
  SDt = (sbeta1 + sgamma1) / 8.0;
  SEt = (sbeta1 - sgamma1) / 8.0;
  beta2m1 = 1.0 - 2.0 * beta;

  ## old code--make sure new code works first */
  ## w2pr2 = (w2 + r2); */
  ## sbetaBA = (rw / 2.0) * (sbeta1 - 2.0 * ((2.0 / t2) - beta1) / (1.0 - 2.0 * beta)); */
  ## transpr[1] = (2.0 * rw / t2) * ((1.0 - R_pow(w2pr2, t1)) / (1 - w2pr2)); */

  s2beta1 = (t1m2 - beta1) / beta2m1;                   ## sum from 1 to t-1 of of (2*beta)^(k-1). */
  transpr[1+1] = rw * s2beta1;                            ## PfB1 = PfDB */
  transpr[6+1] = transpr[1];                              ## PfB0 = PfB1 */

  ## sbetaBA = sum beta1[k] * rw/2 * prob(B->A in remaining steps) */
  sbeta2 = 0.0;
  if(t > 2.0) sbeta2 = (1.0 - beta1 / beta) / (1.0 - beta); ## sum of beta^(k-1) from 1 to (t-2) */
  s2beta2 = (2.0 * t1m2 - (beta1 / beta)) / beta2m1;    ## sum from 1 to t-2 of of (2*beta)^(k-1). */
  sbetaBA = 0.5 * rw * (sbeta2 - s2beta2);

  transpr[0+1] = SDt * w2 + SEt * r2 + sbetaBA;           ## PfA1 = PfDA + PfDEA + PfDBA wrong */
  transpr[5+1] = transpr[0+1];                              ## PfA0 = PfA1 */

  transpr[2+1] = SDt * r2 + SEt * w2 + sbetaBA;           ## PfbC = PfDC + PfDEC + PfDBC wrong*/
  transpr[3+1] = (beta1 + gamma1) / 2.0;                  ## PfD = PfDD */
  transpr[4+1] = (beta1 - gamma1) / 2.0;                  ## PfE = PfDE */

  ## marginal probabilities for one marker */
  transpr[8+1] = -t1 * M_LN2;                             ## Aa */
  transpr[7+1] = log1p(-exp(transpr[8+1])) - M_LN2;         ## AA */
  transpr[9+1] = transpr[7+1];                              ## aa */

  transpr
}
