clear;
close all;

lumFactor = 150;
nBinsF = 30;
nBinsE = 30;

alpha = 0.05;

binMultiplier = 10;
MPos = 30;
MMon = 15;
MCon = 10;

load(['./data/flatData_lumFactor',num2str(lumFactor),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'.mat']);

%% Strict bounds (linear scale)

load(['./results/flatStrictBoundsUnfolded_lumFactor',num2str(lumFactor),'binMultiplier',num2str(binMultiplier),'MPos',num2str(MPos),'MMon',num2str(MMon),'MCon',num2str(MCon),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'alpha',num2str(alpha),'.mat']);

figure;
hold on;

h = bar(0,0,'histc');
set(h,'FaceColor','none','EdgeColor','k');
h = bar(0,0,'histc');
set(h,'FaceColor',[0.90,0.90,1],'EdgeColor',[0.90,0.90,1]);
h = bar(0,0,'histc');
set(h,'FaceColor',[0.75,0.75,1],'EdgeColor',[0.75,0.75,1]);
h = bar(0,0,'histc');
set(h,'FaceColor',[0.45,0.45,1],'EdgeColor',[0.45,0.45,1]);

for i = 1:nBinsE
    rectangle('Position',[binsE(i),lbHPosOc(i)/binWidthsE(i),binWidthsE(i),ubHPosOc(i)/binWidthsE(i)-lbHPosOc(i)/binWidthsE(i)],'FaceColor',[0.90,0.90,1],'EdgeColor','w');
end

for i = 1:nBinsE
    rectangle('Position',[binsE(i),lbHMonOc(i)/binWidthsE(i),binWidthsE(i),ubHMonOc(i)/binWidthsE(i)-lbHMonOc(i)/binWidthsE(i)],'FaceColor',[0.75,0.75,1],'EdgeColor','w');
end

for i = 1:nBinsE
    rectangle('Position',[binsE(i),lbHConOc(i)/binWidthsE(i),binWidthsE(i),ubHConOc(i)/binWidthsE(i)-lbHConOc(i)/binWidthsE(i)],'FaceColor',[0.45,0.45,1],'EdgeColor','w');
end

for i = 1:nBinsE
    line([binsE(i),binsE(i+1)],[fBinsE(i)/binWidthsE(i),fBinsE(i)/binWidthsE(i)],'Color','k','LineStyle','-');
end

plot(gridE,fGridE,'k-');

hold off;

set(gca,'Layer','top');

xlim([lbE,ubE]);
ylim([1400,1900]);
box on;
title('(a) Inclusive jet p_T spectrum, linear scale');
xlabel('Transverse momentum p_T (GeV)');
ylabel('Intensity (1/GeV)');
legend('True','Positive','Positive + decreasing','Positive + decreasing + convex');
set(gcf,'units','centimeters')
set(gcf,'pos',[0 0 15 10])
set(gcf,'paperunits',get(gcf,'units')) 
set(gcf,'paperpos',get(gcf,'pos'))
%print('-depsc2',['./figures/flat_lumFactor',num2str(lumFactor),'binMultiplier',num2str(binMultiplier),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'UnfoldedStrictBoundsLinear.eps']);

%% Compare the lengths of the conservative and unconservative intervals

(ubHPosOc-lbHPosOc)./(ubHPosUc-lbHPosUc)
max((ubHPosOc-lbHPosOc)./(ubHPosUc-lbHPosUc))

(ubHMonOc-lbHMonOc)./(ubHMonUc-lbHMonUc)
max((ubHMonOc-lbHMonOc)./(ubHMonUc-lbHMonUc))

(ubHConOc-lbHConOc)./(ubHConUc-lbHConUc)
max((ubHConOc-lbHConOc)./(ubHConUc-lbHConUc))

%% Coverage studies

nSamples = 1000;

load(['./results/flatStrictBoundsCoverage_lumFactor',num2str(lumFactor),'binMultiplier',num2str(binMultiplier),'nSamples',num2str(nSamples),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'alpha',num2str(alpha),'.mat']);

%disp(coverageBinWisePos);
disp(coverageJointPos);
[~,ci] = binofit(coverageJointPos*nSamples,nSamples);
disp(ci);
%disp(coverageBinWiseMon);
disp(coverageJointMon);
[~,ci] = binofit(coverageJointMon*nSamples,nSamples);
disp(ci);
%disp(coverageBinWiseCon);
disp(coverageJointCon);
[~,ci] = binofit(coverageJointCon*nSamples,nSamples);
disp(ci);