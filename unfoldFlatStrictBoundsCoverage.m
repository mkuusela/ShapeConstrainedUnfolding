clear;
close all;

lumFactor = 150;
nBinsF = 30;
nBinsE = 30;

binMultiplier = 10;
MPos = 30;
MMon = 15;
MCon = 10;

load(['./data/flatData_lumFactor',num2str(lumFactor),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'.mat']);
load(['./results/incJetsStrictBoundsUnfolded_lumFactor',num2str(lumFactor),'binMultiplier',num2str(binMultiplier),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'_init.mat']);

nSamples = 1000;

alpha = 0.05;

lbHPosUc = zeros(nBinsE,nSamples);
lbHPosOc = zeros(nBinsE,nSamples);
lbHMonUc = zeros(nBinsE,nSamples);
lbHMonOc = zeros(nBinsE,nSamples);
lbHConUc = zeros(nBinsE,nSamples);
lbHConOc = zeros(nBinsE,nSamples);
ubHPosUc = zeros(nBinsE,nSamples);
ubHPosOc = zeros(nBinsE,nSamples);
ubHMonUc = zeros(nBinsE,nSamples);
ubHMonOc = zeros(nBinsE,nSamples);
ubHConUc = zeros(nBinsE,nSamples);
ubHConOc = zeros(nBinsE,nSamples);

coverBinWisePos = zeros(nBinsE,nSamples);
coverJointPos = zeros(1,nSamples);
coverBinWiseMon = zeros(nBinsE,nSamples);
coverJointMon = zeros(1,nSamples);
coverBinWiseCon = zeros(nBinsE,nSamples);
coverJointCon = zeros(1,nSamples);

tic;

parfor_progress(nSamples);

parfor iSample = 1:nSamples
    
    [lbHPosUc(:,iSample),lbHPosOc(:,iSample),lbHMonUc(:,iSample),lbHMonOc(:,iSample),lbHConUc(:,iSample),lbHConOc(:,iSample),ubHPosUc(:,iSample),ubHPosOc(:,iSample),ubHMonUc(:,iSample),ubHMonOc(:,iSample),ubHConUc(:,iSample),ubHConOc(:,iSample)] = unfoldStrictBoundsNoConOc(ySeveralSamples(:,iSample),K,KStar,KStarStar,rhoMax,rhoMin,sGrid,Delta,m,binsE,nBinsE,nBinsF,binMultiplier,MPos,MMon,MCon,alpha);
    
    coverBinWisePos(:,iSample) = ubHPosOc(:,iSample) >= fBinsE & lbHPosOc(:,iSample) <= fBinsE;
    coverJointPos(iSample) = all(ubHPosOc(:,iSample) >= fBinsE & lbHPosOc(:,iSample) <= fBinsE);
    coverBinWiseMon(:,iSample) = ubHMonOc(:,iSample) >= fBinsE & lbHMonOc(:,iSample) <= fBinsE;
    coverJointMon(iSample) = all(ubHMonOc(:,iSample) >= fBinsE & lbHMonOc(:,iSample) <= fBinsE);
    % NB: Coverage computed for the unconservatively discretized convex intervals to save CPU time (gives a lower bound for the coverage of the conservatively discretized intervals)
    coverBinWiseCon(:,iSample) = ubHConUc(:,iSample) >= fBinsE & lbHConUc(:,iSample) <= fBinsE;
    coverJointCon(iSample) = all(ubHConUc(:,iSample) >= fBinsE & lbHConUc(:,iSample) <= fBinsE);
    
    parfor_progress;
    
end

parfor_progress(0);

toc;

coverageBinWisePos = mean(coverBinWisePos,2);
coverageJointPos = mean(coverJointPos,2);
coverageBinWiseMon = mean(coverBinWiseMon,2);
coverageJointMon = mean(coverJointMon,2);
coverageBinWiseCon = mean(coverBinWiseCon,2);
coverageJointCon = mean(coverJointCon,2);

save(['./results/flatStrictBoundsCoverage_lumFactor',num2str(lumFactor),'binMultiplier',num2str(binMultiplier),'nSamples',num2str(nSamples),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'alpha',num2str(alpha),'.mat'],'coverageBinWisePos','coverageJointPos','coverageBinWiseMon','coverageJointMon','coverageBinWiseCon','coverageJointCon','lbHPosUc','lbHPosOc','lbHMonUc','lbHMonOc','lbHConUc','lbHConOc','ubHPosUc','ubHPosOc','ubHMonUc','ubHMonOc','ubHConUc','ubHConOc');