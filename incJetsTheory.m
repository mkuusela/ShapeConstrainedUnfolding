% Inclusive jets spectrum in units of pb/GeV
function r = incJetsTheory(pt,L,N0,gamma,alpha,beta,eta)

ECutOff = 3500;

r = L*N0*(pt < ECutOff).*exp(-gamma./pt).*pt.^-alpha.*(1-pt.*cosh(eta)/ECutOff).^beta;