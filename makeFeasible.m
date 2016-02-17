function [nuTildeFeasible,nIter] = makeFeasible(nuTildeStart,cCon,ACon,bCon,nBinsF)

nIter = 0;
nuTildeFeasible = nuTildeStart;
nuTildeFeasible(nuTildeFeasible < 0) = 0; % Clean up negative values
nuTildeFeasible(1:nBinsF) = (nuTildeFeasible(1:nBinsF) > 1e-5).*nuTildeFeasible(1:nBinsF); % Clean up numerically zero positive values

feasible = max(cCon(nuTildeFeasible)) <= 0 && all(ACon*nuTildeFeasible <= bCon);
while ~feasible
    nuTildeFeasible = [0.9999*ones(nBinsF,1); 1.0001*ones(nBinsF,1)].*nuTildeFeasible; % Scale down nuPlus, scale up nuMinus
    feasible = max(cCon(nuTildeFeasible)) <= 0 && all(ACon*nuTildeFeasible <= bCon);
    nIter = nIter + 1;
end