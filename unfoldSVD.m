function [lambdaHat,lambdaHatLbBinWise,lambdaHatUbBinWise,lambdaHatLbJoint,lambdaHatUbJoint] = unfoldSVD(y,K,L,fBinsE,delta,alpha,nBinsE)

alphaPrime = alpha/nBinsE; % Bonferroni correction

CHatInv = diag(1./y);

LTilde = L * diag(1./fBinsE);
KPlus = (K'*CHatInv*K + delta*(LTilde'*LTilde))\(K'*CHatInv);
lambdaHat = KPlus*y;
covLambdaHat = KPlus*diag(y)*KPlus';
lambdaHatUbBinWise = lambdaHat + norminv(1-alpha/2) * sqrt(diag(covLambdaHat));
lambdaHatLbBinWise = lambdaHat - norminv(1-alpha/2) * sqrt(diag(covLambdaHat));
lambdaHatUbJoint = lambdaHat + norminv(1-alphaPrime/2) * sqrt(diag(covLambdaHat));
lambdaHatLbJoint = lambdaHat - norminv(1-alphaPrime/2) * sqrt(diag(covLambdaHat));
