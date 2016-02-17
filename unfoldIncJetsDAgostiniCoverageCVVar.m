clear;
close all;

lumFactor = 150;
nBinsF = 30;
nBinsE = 30;

load(['./data/incJetsData_lumFactor',num2str(lumFactor),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'.mat']);

nSamples = 1000;

alpha = 0.05;

nIterMax = 500; % First try 500 iterations
nIterMaxInc = 20000; % If the minimum of CV not found, increase to 20000 iterations

coverBinWiseDAgostiniTrue = zeros(nBinsE,nSamples);
coverJointDAgostiniTrue = zeros(1,nSamples);

coverBinWiseDAgostiniPert = zeros(nBinsE,nSamples);
coverJointDAgostiniPert = zeros(1,nSamples);

nIterHat = zeros(nSamples,1);

tic;

parfor_progress(nSamples);

parfor iSample = 1:nSamples   

    % D'Agostini (perturbed)
    [~,lambdaHatDAgostiniPertLbBinWise,lambdaHatDAgostiniPertUbBinWise,lambdaHatDAgostiniPertLbJoint,lambdaHatDAgostiniPertUbJoint,nIterHat(iSample)] = unfoldDAgostiniCVVar(ySeveralSamples(:,iSample),KHistPert,fPertBinsE,alpha,nBinsE,nBinsF,nIterMax);
    
    if nIterHat(iSample) == nIterMax % Increase maximum number of iterations
        warning('Retrying with nIterMaxInc');
        [~,lambdaHatDAgostiniPertLbBinWise,lambdaHatDAgostiniPertUbBinWise,lambdaHatDAgostiniPertLbJoint,lambdaHatDAgostiniPertUbJoint,nIterHat(iSample)] = unfoldDAgostiniCVVar(ySeveralSamples(:,iSample),KHistPert,fPertBinsE,alpha,nBinsE,nBinsF,nIterMaxInc);
    end
    
    if nIterHat(iSample) == nIterMaxInc
        warning('Minimum of CV curve not found.');
    end
    
    coverBinWiseDAgostiniPert(:,iSample) = lambdaHatDAgostiniPertUbBinWise >= fBinsE & lambdaHatDAgostiniPertLbBinWise <= fBinsE;
    coverJointDAgostiniPert(iSample) = all(lambdaHatDAgostiniPertUbJoint >= fBinsE & lambdaHatDAgostiniPertLbJoint <= fBinsE);

    % D'Agostini (true)
    [~,lambdaHatDAgostiniTrueLbBinWise,lambdaHatDAgostiniTrueUbBinWise,lambdaHatDAgostiniTrueLbJoint,lambdaHatDAgostiniTrueUbJoint] = unfoldDAgostini(ySeveralSamples(:,iSample),KHistTrue,fBinsE,nIterHat(iSample),alpha,nBinsE,nBinsF);

    coverBinWiseDAgostiniTrue(:,iSample) = lambdaHatDAgostiniTrueUbBinWise >= fBinsE & lambdaHatDAgostiniTrueLbBinWise <= fBinsE;
    coverJointDAgostiniTrue(iSample) = all(lambdaHatDAgostiniTrueUbJoint >= fBinsE & lambdaHatDAgostiniTrueLbJoint <= fBinsE);
    
    parfor_progress;
    
end

parfor_progress(0);

toc;

coverageBinWiseDAgostiniTrue = mean(coverBinWiseDAgostiniTrue,2);
coverageJointDAgostiniTrue = mean(coverJointDAgostiniTrue,2);

coverageBinWiseDAgostiniPert = mean(coverBinWiseDAgostiniPert,2);
coverageJointDAgostiniPert = mean(coverJointDAgostiniPert,2);

save(['./results/incJetsDAgostiniCoverageCVVar_lumFactor',num2str(lumFactor),'nSamples',num2str(nSamples),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'alpha',num2str(alpha),'.mat'],'coverageBinWiseDAgostiniTrue','coverageJointDAgostiniTrue','coverageBinWiseDAgostiniPert','coverageJointDAgostiniPert','coverBinWiseDAgostiniTrue','coverJointDAgostiniTrue','coverBinWiseDAgostiniPert','coverJointDAgostiniPert','nIterHat');