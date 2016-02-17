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

nSamples = 1000;

alpha = 0.05;

coverBinWisePos = zeros(nBinsE,nSamples);
coverJointPos = zeros(1,nSamples);
coverBinWiseMon = zeros(nBinsE,nSamples);
coverJointMon = zeros(1,nSamples);
coverBinWiseCon = zeros(nBinsE,nSamples);
coverJointCon = zeros(1,nSamples);

tic;

%parfor_progress(nSamples);

parfor iSample = 1:nSamples
    
    [lbHPosUc,lbHPosOc,lbHMonUc,lbHMonOc,lbHConUc,lbHConOc,ubHPosUc,ubHPosOc,ubHMonUc,ubHMonOc,ubHConUc,ubHConOc] = unfoldStrictBoundsNoConOc(ySeveralSamples(:,iSample),K,KStar,KStarStar,rhoMax,rhoMin,sGrid,Delta,m,binsE,nBinsE,nBinsF,binMultiplier,MPos,MMon,MCon,alpha);
    
    coverBinWisePos(:,iSample) = ubHPosOc >= fBinsE & lbHPosOc <= fBinsE;
    coverJointPos(iSample) = all(ubHPosOc >= fBinsE & lbHPosOc <= fBinsE);
    coverBinWiseMon(:,iSample) = ubHMonOc >= fBinsE & lbHMonOc <= fBinsE;
    coverJointMon(iSample) = all(ubHMonOc >= fBinsE & lbHMonOc <= fBinsE);
    % NB: Coverage computed for the unconservatively discretized convex intervals to save CPU time (gives a lower bound for the coverage of the conservatively discretized intervals)
    coverBinWiseCon(:,iSample) = ubHConUc >= fBinsE & lbHConUc <= fBinsE;
    coverJointCon(iSample) = all(ubHConUc >= fBinsE & lbHConUc <= fBinsE);
    
    %parfor_progress;
    
end

%parfor_progress(0);

toc;

coverageBinWisePos = mean(coverBinWisePos,2);
coverageJointPos = mean(coverJointPos,2);
coverageBinWiseMon = mean(coverBinWiseMon,2);
coverageJointMon = mean(coverJointMon,2);
coverageBinWiseCon = mean(coverBinWiseCon,2);
coverageJointCon = mean(coverJointCon,2);

save(['./results/incJetsStrictBoundsCoverage_lumFactor',num2str(lumFactor),'binMultiplier',num2str(binMultiplier),'nSamples',num2str(nSamples),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'alpha',num2str(alpha),'.mat'],'coverageBinWisePos','coverageJointPos','coverageBinWiseMon','coverageJointMon','coverageBinWiseCon','coverageJointCon');