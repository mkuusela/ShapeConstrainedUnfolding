clear;
close all;

lumFactor = 150;
nBinsF = 30;
nBinsE = 30;

load(['./data/incJetsData_lumFactor',num2str(lumFactor),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'.mat']);

nSamples = 1000;

alpha = 0.05;

deltas = exp(linspace(log(1e1),log(1e5),20));

L = diag(-2*ones(nBinsE,1)) + diag(ones(nBinsE-1,1),1) + diag(ones(nBinsE-1,1),-1);
L(1,1) = -1;
L(end,end) = -1;

for iDelta=1:length(deltas)

    delta = deltas(iDelta);

    coverBinWiseSVDTrue = zeros(nBinsE,nSamples);
    coverJointSVDTrue = zeros(1,nSamples);

    coverBinWiseSVDPert = zeros(nBinsE,nSamples);
    coverJointSVDPert = zeros(1,nSamples);

    for iSample = 1:nSamples
        
        disp((iDelta-1)*nSamples+iSample);

        % SVD (true)
        [~,lambdaHatSVDTrueLbBinWise,lambdaHatSVDTrueUbBinWise,lambdaHatSVDTrueLbJoint,lambdaHatSVDTrueUbJoint] = unfoldSVD(ySeveralSamples(:,(iDelta-1)*nSamples+iSample),KHistTrue,L,fBinsE,delta,alpha,nBinsE);

        coverBinWiseSVDTrue(:,iSample) = lambdaHatSVDTrueUbBinWise >= fBinsE & lambdaHatSVDTrueLbBinWise <= fBinsE;
        coverJointSVDTrue(iSample) = all(lambdaHatSVDTrueUbJoint >= fBinsE & lambdaHatSVDTrueLbJoint <= fBinsE);

        % SVD (perturbed)
        [~,lambdaHatSVDPertLbBinWise,lambdaHatSVDPertUbBinWise,lambdaHatSVDPertLbJoint,lambdaHatSVDPertUbJoint] = unfoldSVD(ySeveralSamples(:,(iDelta-1)*nSamples+iSample),KHistPert,L,fPertBinsE,delta,alpha,nBinsE);

        coverBinWiseSVDPert(:,iSample) = lambdaHatSVDPertUbBinWise >= fBinsE & lambdaHatSVDPertLbBinWise <= fBinsE;
        coverJointSVDPert(iSample) = all(lambdaHatSVDPertUbJoint >= fBinsE & lambdaHatSVDPertLbJoint <= fBinsE);

    end

    coverageBinWiseSVDTrue = mean(coverBinWiseSVDTrue,2);
    coverageJointSVDTrue = mean(coverJointSVDTrue,2);

    coverageBinWiseSVDPert = mean(coverBinWiseSVDPert,2);
    coverageJointSVDPert = mean(coverJointSVDPert,2);

    save(['./results/incJetsSVDCoverage_lumFactor',num2str(lumFactor),'nSamples',num2str(nSamples),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'alpha',num2str(alpha),'delta',num2str(delta),'.mat'],'coverageBinWiseSVDTrue','coverageJointSVDTrue','coverageBinWiseSVDPert','coverageJointSVDPert');

end