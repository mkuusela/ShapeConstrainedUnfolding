rng(12345);

clear;
close all;

lbE = 400;
ubE = 1000;
nBinsE = 30;
binsE = linspace(lbE,ubE,nBinsE+1)';
binWidthsE = binsE(2:end) - binsE(1:end-1);

lbF = 400;
ubF = 1000;
nBinsF = 30;
binsF = linspace(lbF,ubF,nBinsF+1)';

nGridE = 500;
gridE = linspace(lbE,ubE,nGridE)';
nGridF = 500;
gridF = linspace(lbF,ubF,nGridF)';

% Parameters of the smearing kernel
N = 1;
S = 1;
C = 0.05;

% Parameters of the inclusive jet spectrum
lumFactor = 150; % Size of the 7 TeV dataset
L = lumFactor * 34; % Integrated luminosity in 1/pb
N0 = 1e14;
gamma = 10;
alpha = 5;
beta = 10;
eta = 0;

% Store f
fFun = @(s) incJetsTheory(s,L,N0,gamma,alpha,beta,eta);
fGridE = fFun(gridE);
fBinsE = zeros(nBinsE,1);
for i=1:nBinsE
    fBinsE(i) = quadgk(fFun,binsE(i),binsE(i+1));
end

% Store g = Kf
gGridF = zeros(nGridF,1);
for i=1:nGridF
    t = gridF(i);
    integrand = @(s) normpdf(t,s,s.*sqrt(N^2./s.^2 + S^2./s + C^2)).*fFun(s);
    gGridF(i) = quadgk(integrand,lbE,ubE);
end

% Construct mu
mu = zeros(nBinsF,1);
for i=1:nBinsF
    ki = @(s)(normcdf(binsF(i+1),s,s.*sqrt(N^2./s.^2 + S^2./s + C^2)) - normcdf(binsF(i),s,s.*sqrt(N^2./s.^2 + S^2./s + C^2)));
    fun = @(s)(ki(s).*fFun(s));
    mu(i) = quadgk(fun,lbE,ubE);
end

% Construct K for spline basis (used only for Figure 1)
order = 4;
nKnots = 16;
knots = linspace(lbE,ubE,nKnots+2);
knots = augknt(knots,order);
p = nKnots + order;

KSpline = zeros(nBinsF,p);
for i=1:nBinsF
    ki = @(s)(normcdf(binsF(i+1),s,s.*sqrt(N^2./s.^2 + S^2./s + C^2)) - normcdf(binsF(i),s,s.*sqrt(N^2./s.^2 + S^2./s + C^2)));
    for j=1:p
        beta = zeros(1,p);
        beta(j) = 1;
        Bj = spmak(knots,beta);
        fun = @(s)(ki(s).*fnval(Bj,s));
        intLb = knots(j);
        intUb = knots(j+order);
        KSpline(i,j) = quadgk(fun,intLb,intUb);
    end
end

% Construct K for spline basis without smearing (used to compute the starting point of the spline fit of Figure 1)
KSpline_noSmear = zeros(nBinsF,p);
for i=1:nBinsF
    for j=1:p
        beta = zeros(1,p);
        beta(j) = 1;
        Bj = spmak(knots,beta);
        fun = @(x)fnval(Bj,x);
        intLb = binsF(i);
        intUb = binsF(i+1);
        KSpline_noSmear(i,j) = quadgk(fun,intLb,intUb);
    end
end

% Construct K for histogram discretization using the true spectrum
KHistTrue = zeros(nBinsF,nBinsE);
for i=1:nBinsF
    ki = @(s)(normcdf(binsF(i+1),s,s.*sqrt(N^2./s.^2 + S^2./s + C^2)) - normcdf(binsF(i),s,s.*sqrt(N^2./s.^2 + S^2./s + C^2)));
    for j=1:nBinsE
        fun = @(s)ki(s).*fFun(s);
        KHistTrue(i,j) = quadgk(fun,binsE(j),binsE(j+1))/fBinsE(j);
    end
end

% Construct K for histogram discretization using the perturbed spectrum
N0Pert = 5.5e16;
gammaPert = 10;
alphaPert = 6;
betaPert = 12;

fPertFun = @(s) incJetsTheory(s,L,N0Pert,gammaPert,alphaPert,betaPert,eta);
fPertGridE = fPertFun(gridE);

fPertBinsE = zeros(nBinsE,1);
for i=1:nBinsE
    fPertBinsE(i) = quadgk(fPertFun,binsE(i),binsE(i+1));
end

KHistPert = zeros(nBinsF,nBinsE);
for i=1:nBinsF
    ki = @(s)(normcdf(binsF(i+1),s,s.*sqrt(N^2./s.^2 + S^2./s + C^2)) - normcdf(binsF(i),s,s.*sqrt(N^2./s.^2 + S^2./s + C^2)));
    for j=1:nBinsE
        fun = @(s)ki(s).*fPertFun(s);
        KHistPert(i,j) = quadgk(fun,binsE(j),binsE(j+1))/fPertBinsE(j);
    end
end

% Generate smeared data
y = poissrnd(mu);

% Generate independent replications (used for the coverage studies)
nSamples = 20000;
ySeveralSamples = zeros(nBinsF,nSamples);
for i=1:nSamples
    ySeveralSamples(:,i) = poissrnd(mu);
end

save(['./data/incJetsData_lumFactor',num2str(lumFactor),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'.mat'],'fGridE','fPertGridE','gGridF','fBinsE','fPertBinsE','KSpline','KHistTrue','KHistPert','KSpline_noSmear','mu','y','ySeveralSamples','knots','gridE','gridF','binsE','binWidthsE','binsF','p','nGridE','nGridF','nBinsE','nBinsF','order','lbE','ubE','lbF','ubF');