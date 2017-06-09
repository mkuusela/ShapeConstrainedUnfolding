close all;
clear;

nBinsF = 30;
nBinsE = 30;
binMultiplier = 10;
alpha = 0.05;
nSamples = 1000;

lumFactor = 150;

fileID = fopen(['./results/strictboundsCoverage_alpha',num2str(alpha),'.dat'],'w');

load(['./results/incJetsStrictBoundsCoverage_lumFactor',num2str(lumFactor),'binMultiplier',num2str(binMultiplier),'nSamples',num2str(nSamples),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'alpha',num2str(alpha),'.mat']);

disp(coverageJointPos);
[~,ciJointPos] = binofit(coverageJointPos*nSamples,nSamples);
disp(ciJointPos);
disp(coverageJointMon);
[~,ciJointMon] = binofit(coverageJointMon*nSamples,nSamples);
disp(ciJointMon);
disp(coverageJointCon);
[~,ciJointCon] = binofit(coverageJointCon*nSamples,nSamples);
disp(ciJointCon);

fprintf(fileID,'%s & %1.3f \\textit{(%1.3f, %1.3f)} & %1.3f \\textit{(%1.3f, %1.3f)} & %1.3f \\textit{(%1.3f, %1.3f)} \\\\ \n','Inclusive jets',coverageJointPos,ciJointPos(1),ciJointPos(2),coverageJointMon,ciJointMon(1),ciJointMon(2),coverageJointCon,ciJointCon(1),ciJointCon(2));

load(['./results/linearStrictBoundsCoverage_lumFactor',num2str(lumFactor),'binMultiplier',num2str(binMultiplier),'nSamples',num2str(nSamples),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'alpha',num2str(alpha),'.mat']);

disp(coverageJointPos);
[~,ciJointPos] = binofit(coverageJointPos*nSamples,nSamples);
disp(ciJointPos);
disp(coverageJointMon);
[~,ciJointMon] = binofit(coverageJointMon*nSamples,nSamples);
disp(ciJointMon);
disp(coverageJointCon);
[~,ciJointCon] = binofit(coverageJointCon*nSamples,nSamples);
disp(ciJointCon);

fprintf(fileID,'%s & %1.3f \\textit{(%1.3f, %1.3f)} & %1.3f \\textit{(%1.3f, %1.3f)} & %1.3f \\textit{(%1.3f, %1.3f)} \\\\ \n','Linearly decr.',coverageJointPos,ciJointPos(1),ciJointPos(2),coverageJointMon,ciJointMon(1),ciJointMon(2),coverageJointCon,ciJointCon(1),ciJointCon(2));

load(['./results/flatStrictBoundsCoverage_lumFactor',num2str(lumFactor),'binMultiplier',num2str(binMultiplier),'nSamples',num2str(nSamples),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'alpha',num2str(alpha),'.mat']);

disp(coverageJointPos);
[~,ciJointPos] = binofit(coverageJointPos*nSamples,nSamples);
disp(ciJointPos);
disp(coverageJointMon);
[~,ciJointMon] = binofit(coverageJointMon*nSamples,nSamples);
disp(ciJointMon);
disp(coverageJointCon);
[~,ciJointCon] = binofit(coverageJointCon*nSamples,nSamples);
disp(ciJointCon);

fprintf(fileID,'%s & %1.3f \\textit{(%1.3f, %1.3f)} & %1.3f \\textit{(%1.3f, %1.3f)} & %1.3f \\textit{(%1.3f, %1.3f)} \\\\','Constant',coverageJointPos,ciJointPos(1),ciJointPos(2),coverageJointMon,ciJointMon(1),ciJointMon(2),coverageJointCon,ciJointCon(1),ciJointCon(2));

fclose(fileID);