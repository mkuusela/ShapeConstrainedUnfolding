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

load(['./data/incJetsData_lumFactor',num2str(lumFactor),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'.mat']);

%% Strict bounds (linear scale)

load(['./results/incJetsStrictBoundsUnfolded_lumFactor',num2str(lumFactor),'binMultiplier',num2str(binMultiplier),'MPos',num2str(MPos),'MMon',num2str(MMon),'MCon',num2str(MCon),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'alpha',num2str(alpha),'.mat']);

figure;
hold on;

h = bar(0,0,'style','histc');
set(h,'FaceColor','none','EdgeColor','k');
h = bar(0,0,'style','histc');
set(h,'FaceColor',[0.90,0.90,1],'EdgeColor',[0.90,0.90,1]);
h = bar(0,0,'style','histc');
set(h,'FaceColor',[0.75,0.75,1],'EdgeColor',[0.75,0.75,1]);
h = bar(0,0,'style','histc');
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
ylim([0,1.15*max(fGridE)]);
box on;
title('(a) Inclusive jet p_T spectrum, linear scale');
xlabel('Transverse momentum p_T (GeV)');
ylabel('Intensity (1/GeV)');
legend('True','Positive','Positive + decreasing','Positive + decreasing + convex');
set(gcf,'units','centimeters')
set(gcf,'pos',[0 0 15 10])
set(gcf,'paperunits',get(gcf,'units')) 
set(gcf,'paperpos',get(gcf,'pos'))
print('-depsc2',['./figures/incJets_lumFactor',num2str(lumFactor),'binMultiplier',num2str(binMultiplier),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'UnfoldedStrictBoundsLinear.eps']);

%% Compare the lengths of the conservative and unconservative intervals

(ubHPosOc-lbHPosOc)./(ubHPosUc-lbHPosUc)
max((ubHPosOc-lbHPosOc)./(ubHPosUc-lbHPosUc))

(ubHMonOc-lbHMonOc)./(ubHMonUc-lbHMonUc)
max((ubHMonOc-lbHMonOc)./(ubHMonUc-lbHMonUc))

(ubHConOc-lbHConOc)./(ubHConUc-lbHConUc)
max((ubHConOc-lbHConOc)./(ubHConUc-lbHConUc))

%% Strict bounds (log scale)

load(['./results/incJetsStrictBoundsUnfolded_lumFactor',num2str(lumFactor),'binMultiplier',num2str(binMultiplier),'MPos',num2str(MPos),'MMon',num2str(MMon),'MCon',num2str(MCon),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'alpha',num2str(alpha),'.mat']);

figure;

hold on;

minLb = 1e-5;
lbHPosOc(lbHPosOc <= minLb) = minLb;
lbHMonOc(lbHMonOc <= minLb) = minLb;

h = bar(0,0,'style','histc');
set(h,'FaceColor','none','EdgeColor','k');
h = bar(0,0,'style','histc');
set(h,'FaceColor',[0.90,0.90,1],'EdgeColor',[0.90,0.90,1]);
h = bar(0,0,'style','histc');
set(h,'FaceColor',[0.75,0.75,1],'EdgeColor',[0.75,0.75,1]);
h = bar(0,0,'style','histc');
set(h,'FaceColor',[0.45,0.45,1],'EdgeColor',[0.45,0.45,1]);

for i = 1:nBinsE
    h = rectangle('Position',[binsE(i),lbHPosOc(i)/binWidthsE(i),binWidthsE(i),ubHPosOc(i)/binWidthsE(i)-lbHPosOc(i)/binWidthsE(i)],'FaceColor',[0.90,0.90,1],'EdgeColor','w');
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
ylim([0.7*min(lbHConOc./binWidthsE),1.5*max(ubHPosOc./binWidthsE)]);
box on;
title('(b) Inclusive jet p_T spectrum, log scale');
xlabel('Transverse momentum p_T (GeV)');
ylabel('Intensity (1/GeV)');
legend('True','Positive','Positive + decreasing','Positive + decreasing + convex');
set(gca,'YScale','log');
set(gcf,'units','centimeters')
set(gcf,'pos',[0 0 15 10])
set(gcf,'paperunits',get(gcf,'units')) 
set(gcf,'paperpos',get(gcf,'pos'))
print('-depsc2',['./figures/incJets_lumFactor',num2str(lumFactor),'binMultiplier',num2str(binMultiplier),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'UnfoldedStrictBoundsLog.eps']);

%% SVD, varying delta

nSamples = 1000;

deltas = exp(linspace(log(1e1),log(1e5),20));

jointCoverage = zeros(length(deltas),1);

for iDelta = 1:length(deltas)

    delta = deltas(iDelta);

    load(['./results/incJetsSVDCoverage_lumFactor',num2str(lumFactor),'nSamples',num2str(nSamples),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'alpha',num2str(alpha),'delta',num2str(delta),'.mat']);
    
    jointCoverage(iDelta) = coverageJointSVDPert;

end

[~,ci] = binofit(jointCoverage*nSamples,nSamples);
jointCoverageUb = ci(:,2);
jointCoverageLb = ci(:,1);

figure;
hold on;
line([0.5*deltas(1),2*deltas(end)],[1-alpha,1-alpha],'LineStyle',':','Color','k','LineWidth',1);
errorbar(deltas,jointCoverage,jointCoverage-jointCoverageLb,jointCoverageUb-jointCoverage,'.-b','MarkerSize',7);
hold off;
box on;
title('(a) SVD variant of Tikhonov regularization');
xlabel('Regularization parameter \delta');
ylabel('Simultaneous coverage');
set(gca,'xscale','log');
xlim([1/1.5*deltas(1),1.5*deltas(end)]);
ylim([-0.1,1.1]);
set(gca,'XTick',[1e1,1e2,1e3,1e4,1e5]);
set(gca,'XMinorTick','off');
set(gcf,'units','centimeters')
set(gcf,'pos',[0 0 0.6*15 0.6*10])
set(gcf,'paperunits',get(gcf,'units')) 
set(gcf,'paperpos',get(gcf,'pos'))
print('-depsc2',['./figures/incJets_lumFactor',num2str(lumFactor),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'JointCoverageSVD.eps']);

%% DAgostini, varying nIter

nSamples = 1000;

nIters = [1,3,5,7,9,11,15,20,25,30,40,50];

jointCoverage = zeros(length(nIters),1);

for iIter = 1:length(nIters)

    nIter = nIters(iIter);

    load(['./results/incJetsDAgostiniCoverage_lumFactor',num2str(lumFactor),'nSamples',num2str(nSamples),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'alpha',num2str(alpha),'nIter',num2str(nIter),'.mat']);
    
    jointCoverage(iIter) = coverageJointDAgostiniPert;

end

[~,ci] = binofit(jointCoverage*nSamples,nSamples);
jointCoverageUb = ci(:,2);
jointCoverageLb = ci(:,1);

figure;
hold on;
line([-1,52],[1-alpha,1-alpha],'LineStyle',':','Color','k','LineWidth',1);
errorbar(nIters,jointCoverage,jointCoverage-jointCoverageLb,jointCoverageUb-jointCoverage,'.-b','MarkerSize',7);
hold off;
box on;
title('(b) D''Agostini iteration');
xlabel('Number of iterations');
ylabel('Simultaneous coverage');
xlim([-1,52]);
ylim([-0.1,1.1]);
set(gcf,'units','centimeters')
set(gcf,'pos',[0 0 0.6*15 0.6*10])
set(gcf,'paperunits',get(gcf,'units')) 
set(gcf,'paperpos',get(gcf,'pos'))
print('-depsc2',['./figures/incJets_lumFactor',num2str(lumFactor),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'JointCoverageDAgostini.eps']);


%% SVD, CVVar

nSamples = 1000;

load(['./results/incJetsSVDCoverageCVVar_lumFactor',num2str(lumFactor),'nSamples',num2str(nSamples),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'alpha',num2str(alpha),'.mat']);

[~,ci] = binofit(coverageBinWiseSVDPert*nSamples,nSamples);
binWiseCoverageUb = ci(:,2);
binWiseCoverageLb = ci(:,1);

[~,ci] = binofit(coverageJointSVDPert*nSamples,nSamples);
jointCoverageUb = ci(1,2);
jointCoverageLb = ci(1,1);

figure;
hold on;
line([lbE,ubE],[1-alpha,1-alpha],'LineStyle',':','Color','k','LineWidth',1);
h = bar(binsE(1:end-1),coverageBinWiseSVDPert,'style','histc');
set(h,'FaceColor','none','EdgeColor','b');
errorbar((binsE(1:end-1)+binsE(2:end))/2,coverageBinWiseSVDPert,coverageBinWiseSVDPert-binWiseCoverageLb,binWiseCoverageUb-coverageBinWiseSVDPert,'b','LineStyle','none');
hold off;
box on;
ylim([0.15,1]);
title('(a) SVD, weighted CV');
ylabel('Binwise coverage');
xlabel('Transverse momentum p_T (GeV)');

dim = [0.376 0.25 0.465 0.168];
str = ['Simultaneous coverage: ',num2str(coverageJointSVDPert,'%1.3f'),' (',num2str(jointCoverageLb,'%1.3f'),',',num2str(jointCoverageUb,'%1.3f'),')'];
annotation('textbox',dim,'String',str,'BackgroundColor','w','HorizontalAlignment','center','FontSize',10);

set(gcf,'units','centimeters')
set(gcf,'pos',[0 0 0.6*15 0.6*10])
set(gcf,'paperunits',get(gcf,'units')) 
set(gcf,'paperpos',get(gcf,'pos'))
print('-depsc2',['./figures/incJets_lumFactor',num2str(lumFactor),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'CoverageCVVarSVD.eps']);

%% D'Agostini, CVVar

nSamples = 1000;

load(['./results/incJetsDAgostiniCoverageCVVar_lumFactor',num2str(lumFactor),'nSamples',num2str(nSamples),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'alpha',num2str(alpha),'.mat']);

% Discard runs where the minimum of the CV curve was not found
idx = find(nIterHat == 20000);
coverBinWiseDAgostiniPert(:,idx) = [];
coverJointDAgostiniPert(idx) = [];
coverageBinWiseDAgostiniPert = mean(coverBinWiseDAgostiniPert,2);
coverageJointDAgostiniPert = mean(coverJointDAgostiniPert);
nSamples = nSamples - length(idx);

[~,ci] = binofit(sum(coverBinWiseDAgostiniPert,2),nSamples);
binWiseCoverageUb = ci(:,2);
binWiseCoverageLb = ci(:,1);

[~,ci] = binofit(sum(coverJointDAgostiniPert),nSamples);
jointCoverageUb = ci(1,2);
jointCoverageLb = ci(1,1);

figure;
hold on;
line([lbE,ubE],[1-alpha,1-alpha],'LineStyle',':','Color','k','LineWidth',1);
h = bar(binsE(1:end-1),coverageBinWiseDAgostiniPert,'style','histc');
set(h,'FaceColor','none','EdgeColor','b');
errorbar((binsE(1:end-1)+binsE(2:end))/2,coverageBinWiseDAgostiniPert,coverageBinWiseDAgostiniPert-binWiseCoverageLb,binWiseCoverageUb-coverageBinWiseDAgostiniPert,'b','LineStyle','none');
hold off;
box on;
ylim([0.55,1]);
title('(b) D''Agostini iteration, weighted CV');
ylabel('Binwise coverage');
xlabel('Transverse momentum p_T (GeV)');

dim = [0.376 0.25 0.465 0.168];
str = ['Simultaneous coverage: ',num2str(coverageJointDAgostiniPert,'%1.3f'),' (',num2str(jointCoverageLb,'%1.3f'),',',num2str(jointCoverageUb,'%1.3f'),')'];
annotation('textbox',dim,'String',str,'BackgroundColor','w','HorizontalAlignment','center');

set(gca,'YTick',0.6:0.1:1);

set(gcf,'units','centimeters')
set(gcf,'pos',[0 0 0.6*15 0.6*10])
set(gcf,'paperunits',get(gcf,'units')) 
set(gcf,'paperpos',get(gcf,'pos'))
print('-depsc2',['./figures/incJets_lumFactor',num2str(lumFactor),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'CoverageCVVarDAgostini.eps']);

%% Figure for the introduction

load('./results/splineDemo.mat')

figure;
hold on;

plot(gridE,fGridE,':k','LineWidth',1);
plot(gridE,fHatUncGridE,'-','Color',[0.4,1,0.2],'LineWidth',0.5);
plot(gridE,fHatConGridE,'--','Color',[0.2,0.2,1],'LineWidth',1);

hold off;
ylim([-2500,1.65e4]);
box on;

h = legend('True','Unregularized','Positive + decreasing + convex');
set(h,'FontSize',8);
title('(a) Regularization with shape constraints')

set(gca,'Layer','top');

set(gcf,'units','centimeters')
set(gcf,'pos',[0 0 0.6*15 0.6*10])
set(gcf,'paperunits',get(gcf,'units')) 
set(gcf,'paperpos',get(gcf,'pos'))
print('-depsc2','./figures/incJetsIntro1.eps');


load(['./results/incJetsUnfoldedStrictBounds_lumFactor',num2str(lumFactor),'binMultiplier',num2str(binMultiplier),'MPos',num2str(MPos),'MMon',num2str(MMon),'MCon',num2str(MCon),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'alpha',num2str(alpha),'.mat']);

figure;
hold on;

h = bar(0,0,'style','histc');
set(h,'FaceColor','none','EdgeColor','k');
h = bar(0,0,'style','histc');
set(h,'FaceColor',[0.45,0.45,1],'EdgeColor',[0.45,0.45,1]);

for i = 1:nBinsE
    rectangle('Position',[binsE(i),lbHConOc(i)/binWidthsE(i),binWidthsE(i),ubHConOc(i)/binWidthsE(i)-lbHConOc(i)/binWidthsE(i)],'FaceColor',[0.45,0.45,1],'EdgeColor','none');
end
h = bar(binsE(1:end-1),fBinsE./binWidthsE,'style','histc');
set(h,'FaceColor','none','EdgeColor','k','LineWidth',0.3);
plot(gridE,fGridE,'k');
hold off;
ylim([0,1.5e4]);
box on;

h = legend('True','Positive + decreasing + convex');
set(h,'FontSize',8);
title('(b) Uncertainty quantification with strict bounds')

set(gca,'Layer','top');

set(gcf,'units','centimeters')
set(gcf,'pos',[0 0 0.6*15 0.6*10])
set(gcf,'paperunits',get(gcf,'units')) 
set(gcf,'paperpos',get(gcf,'pos'))
print('-depsc2','./figures/incJetsIntro2.eps');

%% Dual constraints

load(['./results/incJetsStrictBoundsUnfolded_lumFactor',num2str(lumFactor),'binMultiplier',num2str(binMultiplier),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'_init.mat']);
load(['./results/incJetsStrictBoundsUnfolded_lumFactor',num2str(lumFactor),'binMultiplier',num2str(binMultiplier),'MPos',num2str(MPos),'MMon',num2str(MMon),'MCon',num2str(MCon),'nBinsE',num2str(nBinsE),'nBinsF',num2str(nBinsF),'alpha',num2str(alpha),'.mat']);

D = [eye(nBinsF) -eye(nBinsF)];

k = 10;

LPos = double(gridE >= binsE(k) & gridE <= binsE(k+1));
LMon = @(s) max(0,min(s - binsE(k),binsE(k+1)-binsE(k)));
LCon = @(s) 0.5*(s-binsE(k)).^2.*(s > binsE(k)).*(s <= binsE(k+1)) + (0.5*(binsE(k+1)-binsE(k))^2 + (binsE(k+1)-binsE(k))*(s-binsE(k+1))).*(s > binsE(k+1));

optNuLbPosOc = D*optNuTildeLbPosOc(:,k);
optNuUbPosOc = D*optNuTildeUbPosOc(:,k);
optNuLbMonOc = D*optNuTildeLbMonOc(:,k);
optNuUbMonOc = D*optNuTildeUbMonOc(:,k);
optNuLbConOc = D*optNuTildeLbConOc(:,k);
optNuUbConOc = D*optNuTildeUbConOc(:,k);

margin = [0.08,0.1];

figure;
subplot_tight(3,2,1,margin);
plot(gridE,LPos,'k');
hold on;
dashline(sGrid,K*optNuLbPosOc,2,1,2,1,'LineWidth',1.5);
hold off;
box on;
title('(a) Positive, lower bound');
ylim([-0.2,1.2]);
subplot_tight(3,2,2,margin);
plot(gridE,-LPos,'k');
hold on;
dashline(sGrid,K*optNuUbPosOc,2,1,2,1,'LineWidth',1.5);
hold off;
box on;
title('(b) Positive, upper bound');
ylim([-1.2,0.2]);
set(gca,'YTick',[-1,-0.5,0]);
subplot_tight(3,2,3,margin);
plot(sGrid,LMon(sGrid),'k');
hold on;
dashline(sGrid,KStar*optNuLbMonOc,2,1,2,1,'LineWidth',1.5);
hold off;
box on;
title('(c) Decreasing, lower bound');
ylim([-7,25]);
subplot_tight(3,2,4,margin);
plot(sGrid,-LMon(sGrid),'k');
hold on;
dashline(sGrid,KStar*optNuUbMonOc,2,1,2,1,'LineWidth',1.5);
hold off;
box on;
title('(d) Decreasing, upper bound');
ylim([-27,5]);
subplot_tight(3,2,5,margin);
plot(sGrid,LCon(sGrid),'k');
hold on;
dashline(sGrid,KStarStar*optNuLbConOc,2,1,2,1,'LineWidth',1.5);
hold off;
box on;
title('(e) Convex, lower bound');
ylim([-1300,9500]);
h = axes('position', [0.15 0.22 0.15 0.07]);
set(h,'box','on');
plot(sGrid,LCon(sGrid),'k');
xlim([580,590]);
ylim([-10,50]);
hold on;
dashline(sGrid,KStarStar*optNuLbConOc,2,1,2,1,'LineWidth',1.5);
hold off;
set(gca,'FontSize',10);
set(gca,'TickLength',[0.02,0.02]);
subplot_tight(3,2,6,margin);
plot(sGrid,-LCon(sGrid),'k');
hold on;
dashline(sGrid,KStarStar*optNuUbConOc,2,1,2,1,'LineWidth',1.5);
hold off;
box on;
title('(f) Convex, upper bound');
ylim([-9500,1300]);
h = axes('position', [0.62 0.115 0.15 0.07]);
set(h,'box','on');
plot(sGrid,-LCon(sGrid),'k');
xlim([525,645]);
ylim([-950,100]);
hold on;
dashline(sGrid,KStarStar*optNuUbConOc,2,1,2,1,'LineWidth',1.5);
hold off;
set(gca,'XTick',[540,590,640]);
set(gca,'yTick',[-800,-400,0]);
set(gca,'FontSize',10);
set(gca,'TickLength',[0.02,0.02]);

set(gcf,'units','centimeters')
set(gcf,'pos',[0 0 15 18])
set(gcf,'paperunits',get(gcf,'units')) 
set(gcf,'paperpos',get(gcf,'pos'))
print('-depsc2','./figures/constraints.eps');