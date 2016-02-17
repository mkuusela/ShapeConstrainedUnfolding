% Minus log-likelihood of beta (up to a constant that does not depend on beta)
function [ll,grad] = logLikelihood(beta,y,K)

mu = K*beta;
temp = y.*log(mu)-mu;
ll = sum(temp);

grad = K'*(y./mu - 1);

ll = -ll;
grad = -grad;