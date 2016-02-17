function [lambdaHatNew,JacobianNew] = unfoldDAgostiniOneIter(lambdaHatOld,JacobianOld,K,y,nBinsE,nBinsF)

eff = sum(K,1)';

% Form M matrix
denom = K*lambdaHatOld;
M = K./repmat(denom,[1,nBinsE]);

% Calculate next iterate
lambdaHatNew = (lambdaHatOld./eff).*(M'*y);

% Update Jacobian
JacobianNew = zeros(nBinsE,nBinsF);
for j=1:nBinsE
    for m = 1:nBinsF
        v = M(:,j).*y;
        JacobianNew(j,m) = lambdaHatNew(j)/lambdaHatOld(j)*JacobianOld(j,m) + lambdaHatOld(j)/eff(j) * M(m,j) - lambdaHatOld(j)/eff(j)*v'*M*JacobianOld(:,m);
    end
end