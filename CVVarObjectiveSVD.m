function CV = CVVarObjectiveSVD(y,K,L,fBinsE,delta)

CHatInv = diag(1./y);

LTilde = L * diag(1./fBinsE);
KPlus = (K'*CHatInv*K + delta*(LTilde'*LTilde))\(K'*CHatInv);

H = K*KPlus;

muHat = H*y;

CV = norm((1./sqrt(y)).*(y-muHat)./(1-diag(H)))^2;