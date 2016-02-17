clear;
close all;

lumFactor = 150;
nBinsF = 30;
nBinsE = 30;
p = 20;

load(['./data/incJetsData_lumFactor',num2str(lumFactor),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'.mat']);

fun = @(beta) logLikelihood(beta,y,KSpline);

betaInit = computeBetaInit(KSpline_noSmear,y);

grad = 'on';
options = optimoptions('fmincon','Algorithm','sqp','GradObj',grad,'MaxIter',10000,'MaxFunEvals',1000000,'TolX',1e-20,'TolFun',1e-7);

% Unfold without constraints
betaHatUnc = fmincon(fun,betaInit,[],[],[],[],-Inf*ones(p,1),Inf*ones(p,1),[],options);

% Unfold with positivity, monotonicity and convexity constraints
A1 = diag(ones(p,1)) + diag(-ones(p-1,1),1);
A1 = A1(1:end-1,:);
A1 = -A1;
a1 = zeros(p-1,1);

A2 = diag(ones(p,1)) + diag(-2*ones(p-1,1),1) + diag(ones(p-2,1),2);
A2 = A2(1:end-2,:);
A2 = -A2;
a2 = zeros(p-2,1);

betaHatCon = fmincon(fun,betaInit,[A1;A2],[a1;a2],[],[],zeros(p,1),Inf*ones(p,1),[],options);

% Evaluate the splines
fHatUnc = spmak(knots,betaHatUnc');
fHatUncGridE = fnval(fHatUnc,gridE);
fHatCon = spmak(knots,betaHatCon');
fHatConGridE = fnval(fHatCon,gridE);

save('./results/splineDemo.mat','fHatUncGridE','fHatConGridE');