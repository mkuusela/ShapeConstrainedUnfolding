clear;
close all;

lumFactor = 150;
nBinsF = 30;
nBinsE = 30;

load(['./data/incJetsData_lumFactor',num2str(lumFactor),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'.mat']);

% Parameters of the smearing kernel
N = 1;
S = 1;
C = 0.05;

binMultiplier = 10; % Number of sub-bins within each true bin
m = binMultiplier * nBinsE + 1; %NB: this differs by 1 from the definition used in the paper
sGrid = linspace(lbE,ubE,m);
Delta = sGrid(2)-sGrid(1);
K = zeros(m,nBinsF);
KStar = zeros(m,nBinsF);
KStarStar = zeros(m,nBinsF);
rhoMax = zeros(m-1,nBinsF);
rhoMin = zeros(m-1,nBinsF);

for j=1:nBinsF
    
    disp(j);
    
    kj = @(s)(normcdf(binsF(j+1),s,s.*sqrt(N^2./s.^2 + S^2./s + C^2)) - normcdf(binsF(j),s,s.*sqrt(N^2./s.^2 + S^2./s + C^2)));
    kjStar = @(s) quadgk(kj,lbE,s);
    kjStarVectorInput = @(s) arrayfun(kjStar,s);
    kjStarStar = @(s) quadgk(kjStarVectorInput,lbE,s);
    
    tic;
    for i=1:m
        K(i,j) = kj(sGrid(i));   
    end
    toc;
    
    tic;
    for i=1:m
        KStar(i,j) = kjStar(sGrid(i));
    end
    toc;
    
    tic;
    for i=1:m
        KStarStar(i,j) = kjStarStar(sGrid(i));
    end
    toc;
    
    tic;
    for i=1:m-1
        fun = @(s) -kj(s);
        [~,val] = fminbnd(fun,sGrid(i),sGrid(i+1));
        rhoMax(i,j) = -val;
    end
    toc;
    
    tic;
    for i=1:m-1
        fun = @(s) kj(s);
        [~,val] = fminbnd(fun,sGrid(i),sGrid(i+1));
        rhoMin(i,j) = val;
    end
    toc;

end

save(['./results/incJetsStrictBoundsUnfolded_lumFactor',num2str(lumFactor),'binMultiplier',num2str(binMultiplier),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'_init.mat'],'K','KStar','KStarStar','rhoMax','rhoMin','sGrid','Delta','m');