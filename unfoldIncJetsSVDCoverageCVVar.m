clear;
close all;

lumFactor = 150;
nBinsF = 30;
nBinsE = 30;

load(['./data/incJetsData_lumFactor',num2str(lumFactor),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'.mat']);

nSamples = 1000;

alpha = 0.05;

L = diag(-2*ones(nBinsE,1)) + diag(ones(nBinsE-1,1),1) + diag(ones(nBinsE-1,1),-1);
L(1,1) = -1;
L(end,end) = -1;

coverBinWiseSVDTrue = zeros(nBinsE,nSamples);
coverJointSVDTrue = zeros(1,nSamples);

coverBinWiseSVDPert = zeros(nBinsE,nSamples);
coverJointSVDPert = zeros(1,nSamples);

deltaHat = zeros(nSamples,1);

for iSample = 1:nSamples

    % Optimize delta
    fun = @(delta) CVVarObjectiveSVD(ySeveralSamples(:,iSample),KHistPert,L,fPertBinsE,delta);
    deltaHat(iSample) = fminbnd(fun,1e-3,1e7);
    disp(deltaHat(iSample));
    
    % SVD (true)
    [~,lambdaHatSVDTrueLbBinWise,lambdaHatSVDTrueUbBinWise,lambdaHatSVDTrueLbJoint,lambdaHatSVDTrueUbJoint] = unfoldSVD(ySeveralSamples(:,iSample),KHistTrue,L,fBinsE,deltaHat(iSample),alpha,nBinsE);

    coverBinWiseSVDTrue(:,iSample) = lambdaHatSVDTrueUbBinWise >= fBinsE & lambdaHatSVDTrueLbBinWise <= fBinsE;
    coverJointSVDTrue(iSample) = all(lambdaHatSVDTrueUbJoint >= fBinsE & lambdaHatSVDTrueLbJoint <= fBinsE);
    
    % SVD (perturbed)
    [~,lambdaHatSVDPertLbBinWise,lambdaHatSVDPertUbBinWise,lambdaHatSVDPertLbJoint,lambdaHatSVDPertUbJoint] = unfoldSVD(ySeveralSamples(:,iSample),KHistPert,L,fPertBinsE,deltaHat(iSample),alpha,nBinsE);

    coverBinWiseSVDPert(:,iSample) = lambdaHatSVDPertUbBinWise >= fBinsE & lambdaHatSVDPertLbBinWise <= fBinsE;
    coverJointSVDPert(iSample) = all(lambdaHatSVDPertUbJoint >= fBinsE & lambdaHatSVDPertLbJoint <= fBinsE);
    
end

coverageBinWiseSVDTrue = mean(coverBinWiseSVDTrue,2);
coverageJointSVDTrue = mean(coverJointSVDTrue,2);

coverageBinWiseSVDPert = mean(coverBinWiseSVDPert,2);
coverageJointSVDPert = mean(coverJointSVDPert,2);

save(['./results/incJetsSVDCoverageCVVar_lumFactor',num2str(lumFactor),'nSamples',num2str(nSamples),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'alpha',num2str(alpha),'.mat'],'coverageBinWiseSVDTrue','coverageJointSVDTrue','coverageBinWiseSVDPert','coverageJointSVDPert','deltaHat');