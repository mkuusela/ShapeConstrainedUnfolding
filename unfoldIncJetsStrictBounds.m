clear;
close all;

lumFactor = 150;
nBinsF = 30;
nBinsE = 30;

binMultiplier = 10;
MPos = 30;
MMon = 15;
MCon = 10;

load(['./data/incJetsData_lumFactor',num2str(lumFactor),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'.mat']);
load(['./results/incJetsStrictBoundsUnfolded_lumFactor',num2str(lumFactor),'binMultiplier',num2str(binMultiplier),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'_init.mat']);

alpha = 0.05;

tic;
[lbHPosUc,lbHPosOc,lbHMonUc,lbHMonOc,lbHConUc,lbHConOc,ubHPosUc,ubHPosOc,ubHMonUc,ubHMonOc,ubHConUc,ubHConOc,optNuTildeLbPosUc,optNuTildeLbPosOc,optNuTildeLbMonUc,optNuTildeLbMonOc,optNuTildeLbConUc,optNuTildeLbConOc,optNuTildeUbPosUc,optNuTildeUbPosOc,optNuTildeUbMonUc,optNuTildeUbMonOc,optNuTildeUbConUc,optNuTildeUbConOc] = unfoldStrictBounds(y,K,KStar,KStarStar,rhoMax,rhoMin,sGrid,Delta,m,binsE,nBinsE,nBinsF,binMultiplier,MPos,MMon,MCon,alpha);
toc;

save(['./results/incJetsStrictBoundsUnfolded_lumFactor',num2str(lumFactor),'binMultiplier',num2str(binMultiplier),'MPos',num2str(MPos),'MMon',num2str(MMon),'MCon',num2str(MCon),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'alpha',num2str(alpha),'.mat'],'lbHPosUc','lbHPosOc','lbHMonUc','lbHMonOc','lbHConUc','lbHConOc','ubHPosUc','ubHPosOc','ubHMonUc','ubHMonOc','ubHConUc','ubHConOc','optNuTildeLbPosUc','optNuTildeLbPosOc','optNuTildeLbMonUc','optNuTildeLbMonOc','optNuTildeLbConUc','optNuTildeLbConOc','optNuTildeUbPosUc','optNuTildeUbPosOc','optNuTildeUbMonUc','optNuTildeUbMonOc','optNuTildeUbConUc','optNuTildeUbConOc');

%% Sanity checks

disp(lbHPosOc./lbHPosUc);
disp(ubHPosOc./ubHPosUc);

disp(lbHMonOc./lbHMonUc);
disp(ubHMonOc./ubHMonUc);

disp(lbHConOc./lbHConUc);
disp(ubHConOc./ubHConUc);

disp(max(max(abs(optNuTildeLbPosUc))));
disp(max(max(abs(optNuTildeUbPosUc))));
disp(max(max(abs(optNuTildeLbPosOc))));
disp(max(max(abs(optNuTildeUbPosOc))));

disp(max(max(abs(optNuTildeLbMonUc))));
disp(max(max(abs(optNuTildeUbMonUc))));
disp(max(max(abs(optNuTildeLbMonOc))));
disp(max(max(abs(optNuTildeUbMonOc))));

disp(max(max(abs(optNuTildeLbConUc))));
disp(max(max(abs(optNuTildeUbConUc))));
disp(max(max(abs(optNuTildeLbConOc))));
disp(max(max(abs(optNuTildeUbConOc))));