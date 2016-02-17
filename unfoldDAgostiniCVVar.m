function [lambdaHat,lambdaHatLbBinWise,lambdaHatUbBinWise,lambdaHatLbJoint,lambdaHatUbJoint,nIterHat] = unfoldDAgostiniCVVar(y,K,fBinsE,alpha,nBinsE,nBinsF,nIterMax)

alphaPrime = alpha/nBinsE; % Bonferroni correction

predError = zeros(nBinsF,nIterMax);

for leaveOutIdx = 1:nBinsF
    
    KLeaveOut = K;
    KLeaveOut(leaveOutIdx,:) = [];
    yLeaveOut = y;
    yLeaveOut(leaveOutIdx) = [];
    
    lambdaHat = fBinsE;
    Jacobian = zeros(nBinsE,nBinsF-1);
    
    for iIter = 1:nIterMax
        [lambdaHat,Jacobian] = unfoldDAgostiniOneIter(lambdaHat,Jacobian,KLeaveOut,yLeaveOut,nBinsE,nBinsF-1);
        
        predError(leaveOutIdx,iIter) = (K(leaveOutIdx,:)*lambdaHat-y(leaveOutIdx))^2;
    end
    
end

predError = repmat((1./y),[1,nIterMax]).*predError;

CV = sum(predError,1);

[~,nIterHat] = min(CV);

if nIterHat == nIterMax
    warning('Maximum number of iterations reached.')
end

lambdaHat = fBinsE;
Jacobian = zeros(nBinsE,nBinsF);

for iIter = 1:nIterHat
    [lambdaHat,Jacobian] = unfoldDAgostiniOneIter(lambdaHat,Jacobian,K,y,nBinsE,nBinsF);
end

covLambdaHat = Jacobian*diag(y)*Jacobian';
lambdaHatUbBinWise = lambdaHat + norminv(1-alpha/2) * sqrt(diag(covLambdaHat));
lambdaHatLbBinWise = lambdaHat - norminv(1-alpha/2) * sqrt(diag(covLambdaHat));
lambdaHatUbJoint = lambdaHat + norminv(1-alphaPrime/2) * sqrt(diag(covLambdaHat));
lambdaHatLbJoint = lambdaHat - norminv(1-alphaPrime/2) * sqrt(diag(covLambdaHat));
