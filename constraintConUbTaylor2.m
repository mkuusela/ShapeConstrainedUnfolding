function [r,req] = constraintConUbTaylor2(nuTilde,KStar,KStarStar,rhoMax,rhoMin,D,sGrid,Ekmin,Ekmax,m)

r1 = zeros(m-1,1);
r2 = zeros(m-1,1);
r3 = zeros(m-1,1);

for i=1:m-1
    
    A = KStarStar(i,:)*D*nuTilde;
    B = KStar(i,:)*D*nuTilde;
    C = 0.5*[rhoMax(i,:) -rhoMin(1,:)]*nuTilde;
    
    if sGrid(i) < Ekmin % Flat part
        ATilde = -C;
        BTilde = 2*C*sGrid(i)-B;
        CTilde = -A+B*sGrid(i)-C*sGrid(i)^2;
    elseif sGrid(i) < Ekmax % Quadratic part
        ATilde = -0.5-C;
        BTilde = 2*C*sGrid(i)+Ekmin-B;
        CTilde = -0.5*Ekmin^2-A+B*sGrid(i)-C*sGrid(i)^2;
    else % Linear part
        ATilde = -C;
        BTilde = -Ekmax+Ekmin+2*C*sGrid(i)-B;
        CTilde = -0.5*(Ekmax-Ekmin)^2+Ekmax*(Ekmax-Ekmin)-A+B*sGrid(i)-C*sGrid(i)^2;
    end

    sStar = -BTilde/(2*ATilde);

    r1(i) = -(ATilde*sGrid(i)^2 + BTilde*sGrid(i) + CTilde);
    r2(i) = -(ATilde*sGrid(i+1)^2 + BTilde*sGrid(i+1) + CTilde);
    if sStar > sGrid(i) && sStar < sGrid(i+1)
        r3(i) = -(ATilde*sStar^2 + BTilde*sStar + CTilde);
    else
        r3(i) = 0;
    end
    
end

%r = [r1; r2; r3];
r = [r2; r3];
req = [];