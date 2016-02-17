% Shape constrained strict bounds confidence intervals
% Uc = potential undercoverage (unconservative discretization)
% Oc = potential overcoverage (conservative discretization)

function [lbHPosUc,lbHPosOc,lbHMonUc,lbHMonOc,lbHConUc,lbHConOc,ubHPosUc,ubHPosOc,ubHMonUc,ubHMonOc,ubHConUc,ubHConOc,optNuTildeLbPosUc,optNuTildeLbPosOc,optNuTildeLbMonUc,optNuTildeLbMonOc,optNuTildeLbConUc,optNuTildeLbConOc,optNuTildeUbPosUc,optNuTildeUbPosOc,optNuTildeUbMonUc,optNuTildeUbMonOc,optNuTildeUbConUc,optNuTildeUbConOc] = unfoldStrictBoundsNoConOc(y,K,KStar,KStarStar,rhoMax,rhoMin,sGrid,Delta,m,binsE,nBinsE,nBinsF,binMultiplier,MPos,MMon,MCon,alpha)

alphaPrime = 1-(1-alpha)^(1/nBinsF);

muLb = 0.5 * chi2inv(alphaPrime/2,2*y);
muUb = 0.5 * chi2inv(1-alphaPrime/2,2*(y+1));

yTilde = (muLb + muUb)/2;
l = muUb - yTilde;

lbHPosUc = zeros(nBinsE,1);
lbHPosOc = zeros(nBinsE,1);
lbHMonUc = zeros(nBinsE,1);
lbHMonOc = zeros(nBinsE,1);
lbHConUc = zeros(nBinsE,1);
lbHConOc = zeros(nBinsE,1);

ubHPosUc = zeros(nBinsE,1);
ubHPosOc = zeros(nBinsE,1);
ubHMonUc = zeros(nBinsE,1);
ubHMonOc = zeros(nBinsE,1);
ubHConUc = zeros(nBinsE,1);
ubHConOc = zeros(nBinsE,1);

optNuTildeLbPosUc = zeros(2*nBinsF,nBinsE);
optNuTildeLbPosOc = zeros(2*nBinsF,nBinsE);
optNuTildeLbMonUc = zeros(2*nBinsF,nBinsE);
optNuTildeLbMonOc = zeros(2*nBinsF,nBinsE);
optNuTildeLbConUc = zeros(2*nBinsF,nBinsE);
optNuTildeLbConOc = zeros(2*nBinsF,nBinsE);

optNuTildeUbPosUc = zeros(2*nBinsF,nBinsE);
optNuTildeUbPosOc = zeros(2*nBinsF,nBinsE);
optNuTildeUbMonUc = zeros(2*nBinsF,nBinsE);
optNuTildeUbMonOc = zeros(2*nBinsF,nBinsE);
optNuTildeUbConUc = zeros(2*nBinsF,nBinsE);
optNuTildeUbConOc = zeros(2*nBinsF,nBinsE);

for k=1:nBinsE
    
    % Compute the strict bounds for the kth true bin
    
    disp(' ');
    disp('-----------------------------');
    disp(['Bin number: ',num2str(k)]);
    disp('-----------------------------');
    
    D = [eye(nBinsF) -eye(nBinsF)];
    lTilde = [l; l];
    f = lTilde - D'*yTilde;
    
    fCon = @(nuTilde) funCon(nuTilde,f);
    zeroFun = @(nuTilde) zeros(2*m-2,1);
    
    % Unconservative discretization
    APosUc = K*D;
    if k ~= nBinsE
        bPosUc = double(sGrid >= binsE(k) & sGrid < binsE(k+1))';
    else
        bPosUc = double(sGrid >= binsE(k) & sGrid <= binsE(k+1))';
    end

    AMonUc = KStar*D;
    LMon = @(s) max(0,min(s - binsE(k),binsE(k+1)-binsE(k)));
    bMonUc = LMon(sGrid');
    
    AConUc = [KStarStar*D; KStar(end,:)*D];
    LCon = @(s) 0.5*(s-binsE(k)).^2.*(s > binsE(k)).*(s <= binsE(k+1)) + (0.5*(binsE(k+1)-binsE(k))^2 + (binsE(k+1)-binsE(k))*(s-binsE(k+1))).*(s > binsE(k+1));
    bConUc = [LCon(sGrid'); binsE(k+1)-binsE(k)];

    % Conservative discretization
    APosOc = [rhoMax -rhoMin];
    bPosOc = bPosUc(1:m-1);
    
    AMonOc = KStar(1:m-1,:)*D + [rhoMax*Delta -rhoMin*Delta];
    bMonOcLb = LMon(sGrid(2:m)');
    bMonOcUb = -LMon(sGrid(2:m)');
    
    cConLb = @(nuTilde) constraintConLbTaylor2(nuTilde,KStar,KStarStar,rhoMax,rhoMin,D,sGrid,binsE(k),binsE(k+1),m);
    cConUb = @(nuTilde) constraintConUbTaylor2(nuTilde,KStar,KStarStar,rhoMax,rhoMin,D,sGrid,binsE(k),binsE(k+1),m);
    AConOc = KStar(end,:)*D;
    bConOcLb = binsE(k+1)-binsE(k);
    bConOcUb = -(binsE(k+1)-binsE(k));
    
    optionsLinprog = optimoptions('linprog','Algorithm','interior-point','MaxIter',1000,'TolFun',1e-6);
    optionsFmincon = optimoptions('fmincon','MaxFunEvals',10000,'GradObj','off','Algorithm','sqp','Display','iter-detailed','TolFun',1e-3);

    % Compute the lower bounds    
    [optNuTildeLbPosUc(:,k), lbHPosUc(k), exitflag, output] = linprog(f,APosUc,bPosUc,[],[],zeros(2*nBinsF,1),MPos*ones(2*nBinsF,1),[],optionsLinprog);
    lbHPosUc(k) = -lbHPosUc(k);
    disp(lbHPosUc(k));

    [optNuTildeLbPosOcTemp, ~, exitflag, output] = linprog(f,APosOc,bPosOc,[],[],zeros(2*nBinsF,1),MPos*ones(2*nBinsF,1),[],optionsLinprog);
    disp(['Bound before scaling: ',num2str(-f'*optNuTildeLbPosOcTemp)]);
    disp(['Maximum constraint violation: ',num2str(max(APosOc*optNuTildeLbPosOcTemp - bPosOc))]);
    [optNuTildeLbPosOcTemp,nIter] = makeFeasible(optNuTildeLbPosOcTemp,zeroFun,APosOc,bPosOc,nBinsF);
    lbHPosOc(k) = -f'*optNuTildeLbPosOcTemp;
    optNuTildeLbPosOc(:,k) = optNuTildeLbPosOcTemp;
    disp(lbHPosOc(k));
    disp(['Number of scaling iterations: ',num2str(nIter)]);

    [optNuTildeLbMonUc(:,k), lbHMonUc(k), exitflag, output] = linprog(f,AMonUc,bMonUc,[],[],zeros(2*nBinsF,1),MMon*ones(2*nBinsF,1),[],optionsLinprog);
    lbHMonUc(k) = -lbHMonUc(k);
    disp(lbHMonUc(k));

    [optNuTildeLbMonOcTemp, ~, exitflag, output] = linprog(f,AMonOc,bMonOcLb,[],[],zeros(2*nBinsF,1),MMon*ones(2*nBinsF,1),[],optionsLinprog);
    disp(['Bound before scaling: ',num2str(-f'*optNuTildeLbMonOcTemp)]);
    disp(['Maximum constraint violation: ',num2str(max(AMonOc*optNuTildeLbMonOcTemp - bMonOcLb))]);
    [optNuTildeLbMonOcTemp,nIter] = makeFeasible(optNuTildeLbMonOcTemp,zeroFun,AMonOc,bMonOcLb,nBinsF);
    lbHMonOc(k) = -f'*optNuTildeLbMonOcTemp;
    optNuTildeLbMonOc(:,k) = optNuTildeLbMonOcTemp;
    disp(lbHMonOc(k));
    disp(['Number of scaling iterations: ',num2str(nIter)]);
 
    [optNuTildeLbConUc(:,k), lbHConUc(k), exitflag, output] = linprog(f,AConUc,bConUc,[],[],zeros(2*nBinsF,1),MCon*ones(2*nBinsF,1),[],optionsLinprog);
    lbHConUc(k) = -lbHConUc(k);
    disp(lbHConUc(k));
        
%     [initNuTildeLbConOc,nIter] = makeFeasible(optNuTildeLbConUc(:,k),cConLb,AConOc,bConOcLb,nBinsF);
%     disp(['Bound after scaling: ',num2str(-f'*initNuTildeLbConOc)]);
%     disp(['Number of scaling iterations: ',num2str(nIter)]);
%     [optNuTildeLbConOcTemp, ~, exitflag, output] = fmincon(fCon,initNuTildeLbConOc,AConOc,bConOcLb,[],[],zeros(2*nBinsF,1),MCon*ones(2*nBinsF,1),cConLb,optionsFmincon);
%     [optNuTildeLbConOcTemp,nIter] = makeFeasible(optNuTildeLbConOcTemp,cConLb,AConOc,bConOcLb,nBinsF);
%     if -f'*initNuTildeLbConOc <= -f'*optNuTildeLbConOcTemp
%         lbHConOc(k) = -f'*optNuTildeLbConOcTemp;
%         optNuTildeLbConOc(:,k) = optNuTildeLbConOcTemp;
%     else % SQP iteration worsened the bound
%         warning('Reverting to the starting point.')
%         lbHConOc(k) = -f'*initNuTildeLbConOc;
%         optNuTildeLbConOc(:,k) = initNuTildeLbConOc;
%     end
%     disp(lbHConOc(k));
%     disp(['Number of scaling iterations: ',num2str(nIter)]);
    
    % Compute the upper bounds
    [optNuTildeUbPosUc(:,k), ubHPosUc(k), exitflag, output] = linprog(f,APosUc,-bPosUc,[],[],zeros(2*nBinsF,1),MPos*ones(2*nBinsF,1),[],optionsLinprog);
    disp(ubHPosUc(k));
    
    [optNuTildeUbPosOcTemp, ~, exitflag, output] = linprog(f,APosOc,-bPosOc,[],[],zeros(2*nBinsF,1),MPos*ones(2*nBinsF,1),[],optionsLinprog);
    disp(['Bound before scaling: ',num2str(f'*optNuTildeUbPosOcTemp)]);
    disp(['Maximum constraint violation: ',num2str(max(APosOc*optNuTildeUbPosOcTemp + bPosOc))]);
    [optNuTildeUbPosOcTemp,nIter] = makeFeasible(optNuTildeUbPosOcTemp,zeroFun,APosOc,-bPosOc,nBinsF);
    ubHPosOc(k) = f'*optNuTildeUbPosOcTemp;
    optNuTildeUbPosOc(:,k) = optNuTildeUbPosOcTemp;
    disp(ubHPosOc(k));
    disp(['Number of scaling iterations: ',num2str(nIter)]);

    [optNuTildeUbMonUc(:,k), ubHMonUc(k), exitflag, output] = linprog(f,AMonUc,-bMonUc,[],[],zeros(2*nBinsF,1),MMon*ones(2*nBinsF,1),[],optionsLinprog);
    disp(ubHMonUc(k));

    [optNuTildeUbMonOcTemp, ~, exitflag, output] = linprog(f,AMonOc,bMonOcUb,[],[],zeros(2*nBinsF,1),MMon*ones(2*nBinsF,1),[],optionsLinprog);
    disp(['Bound before scaling: ',num2str(f'*optNuTildeUbMonOcTemp)]);
    disp(['Maximum constraint violation: ',num2str(max(AMonOc*optNuTildeUbMonOcTemp - bMonOcUb))]);
    [optNuTildeUbMonOcTemp,nIter] = makeFeasible(optNuTildeUbMonOcTemp,zeroFun,AMonOc,bMonOcUb,nBinsF);
    ubHMonOc(k) = f'*optNuTildeUbMonOcTemp;
    optNuTildeUbMonOc(:,k) = optNuTildeUbMonOcTemp;
    disp(ubHMonOc(k));
    disp(['Number of scaling iterations: ',num2str(nIter)]);

    [optNuTildeUbConUc(:,k), ubHConUc(k), exitflag, output] = linprog(f,AConUc,-bConUc,[],[],zeros(2*nBinsF,1),MCon*ones(2*nBinsF,1),[],optionsLinprog);
    disp(ubHConUc(k));

%     [initNuTildeUbConOc,nIter] = makeFeasible(optNuTildeUbConUc(:,k),cConUb,AConOc,bConOcUb,nBinsF);
%     disp(['Bound after scaling: ',num2str(f'*initNuTildeUbConOc)]);
%     disp(['Number of scaling iterations: ',num2str(nIter)]);
%     [optNuTildeUbConOcTemp, ~, exitflag, output] = fmincon(fCon,initNuTildeUbConOc,AConOc,bConOcUb,[],[],zeros(2*nBinsF,1),MCon*ones(2*nBinsF,1),cConUb,optionsFmincon);
%     [optNuTildeUbConOcTemp,nIter] = makeFeasible(optNuTildeUbConOcTemp,cConUb,AConOc,bConOcUb,nBinsF);
%     if f'*initNuTildeUbConOc >= f'*optNuTildeUbConOcTemp
%         ubHConOc(k) = f'*optNuTildeUbConOcTemp;
%         optNuTildeUbConOc(:,k) = optNuTildeUbConOcTemp;
%     else % SQP iteration worsened the bound
%         warning('Reverting to the starting point.')
%         ubHConOc(k) = f'*initNuTildeUbConOc;
%         optNuTildeUbConOc(:,k) = initNuTildeUbConOc;
%     end
%     disp(ubHConOc(k));
%     disp(['Number of scaling iterations: ',num2str(nIter)]);

end