function [lambdaHat,lambdaHatLbBinWise,lambdaHatUbBinWise,lambdaHatLbJoint,lambdaHatUbJoint] = unfoldDAgostini(y,K,fBinsE,nIter,alpha,nBinsE,nBinsF)

alphaPrime = alpha/nBinsE; % Bonferroni correction

lambdaHat = fBinsE; % Starting point
Jacobian = zeros(nBinsE,nBinsF);

for iIter = 1:nIter
    [lambdaHat,Jacobian] = unfoldDAgostiniOneIter(lambdaHat,Jacobian,K,y,nBinsE,nBinsF);
end

covLambdaHat = Jacobian*diag(y)*Jacobian';
lambdaHatUbBinWise = lambdaHat + norminv(1-alpha/2) * sqrt(diag(covLambdaHat));
lambdaHatLbBinWise = lambdaHat - norminv(1-alpha/2) * sqrt(diag(covLambdaHat));
lambdaHatUbJoint = lambdaHat + norminv(1-alphaPrime/2) * sqrt(diag(covLambdaHat));
lambdaHatLbJoint = lambdaHat - norminv(1-alphaPrime/2) * sqrt(diag(covLambdaHat));
