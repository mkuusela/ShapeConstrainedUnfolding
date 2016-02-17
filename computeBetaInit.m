function betaInit = computeBetaInit(K_noSmear,y)

betaInit = lsqnonneg(K_noSmear,y);
minPos = min(betaInit(betaInit ~= 0));
nZeros = sum(betaInit == 0);
betaInit(betaInit == 0) = minPos * ones(nZeros,1);