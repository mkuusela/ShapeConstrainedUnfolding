clear;
close all;

lumFactor = 150;
nBinsF = 30;
nBinsE = 30;

load(['./data/incJetsData_lumFactor',num2str(lumFactor),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'.mat']);

nSamples = 1000;

alpha = 0.05;

nIters = [1,3,5,7,9,11,15,20,25,30,40,50];

for iIter = 1:length(nIters)

    nIter = nIters(iIter);
    
    disp(nIter);

    coverBinWiseDAgostiniTrue = zeros(nBinsE,nSamples);
    coverJointDAgostiniTrue = zeros(1,nSamples);

    coverBinWiseDAgostiniPert = zeros(nBinsE,nSamples);
    coverJointDAgostiniPert = zeros(1,nSamples);
    
    tic;

    parfor iSample = 1:nSamples

        % D'Agostini (true)
        [~,lambdaHatDAgostiniTrueLbBinWise,lambdaHatDAgostiniTrueUbBinWise,lambdaHatDAgostiniTrueLbJoint,lambdaHatDAgostiniTrueUbJoint] = unfoldDAgostini(ySeveralSamples(:,(iIter-1)*nSamples+iSample),KHistTrue,fBinsE,nIter,alpha,nBinsE,nBinsF);

        coverBinWiseDAgostiniTrue(:,iSample) = lambdaHatDAgostiniTrueUbBinWise >= fBinsE & lambdaHatDAgostiniTrueLbBinWise <= fBinsE;
        coverJointDAgostiniTrue(iSample) = all(lambdaHatDAgostiniTrueUbJoint >= fBinsE & lambdaHatDAgostiniTrueLbJoint <= fBinsE);

        % D'Agostini (perturbed)
        [~,lambdaHatDAgostiniPertLbBinWise,lambdaHatDAgostiniPertUbBinWise,lambdaHatDAgostiniPertLbJoint,lambdaHatDAgostiniPertUbJoint] = unfoldDAgostini(ySeveralSamples(:,(iIter-1)*nSamples+iSample),KHistPert,fPertBinsE,nIter,alpha,nBinsE,nBinsF);

        coverBinWiseDAgostiniPert(:,iSample) = lambdaHatDAgostiniPertUbBinWise >= fBinsE & lambdaHatDAgostiniPertLbBinWise <= fBinsE;
        coverJointDAgostiniPert(iSample) = all(lambdaHatDAgostiniPertUbJoint >= fBinsE & lambdaHatDAgostiniPertLbJoint <= fBinsE);

    end
    
    toc;

    coverageBinWiseDAgostiniTrue = mean(coverBinWiseDAgostiniTrue,2);
    coverageJointDAgostiniTrue = mean(coverJointDAgostiniTrue,2);

    coverageBinWiseDAgostiniPert = mean(coverBinWiseDAgostiniPert,2);
    coverageJointDAgostiniPert = mean(coverJointDAgostiniPert,2);

    save(['./results/incJetsDAgostiniCoverage_lumFactor',num2str(lumFactor),'nSamples',num2str(nSamples),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'alpha',num2str(alpha),'nIter',num2str(nIter),'.mat'],'coverageBinWiseDAgostiniTrue','coverageJointDAgostiniTrue','coverageBinWiseDAgostiniPert','coverageJointDAgostiniPert');

end