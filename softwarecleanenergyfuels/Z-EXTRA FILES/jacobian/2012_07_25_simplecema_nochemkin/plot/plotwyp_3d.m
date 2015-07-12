clc;clear;
chem=ckinit;
P=1:5:100; %atm then converted below to erg/cm^3
T=300:500:3000;
X = zeros(9,1);
X =X+0.01;
Phi=1;

%initial mole fractions and mass fractions
X(1) = Phi*2;
X(2) = 1;
X(9) = 3.76;
X = X/sum(X);
W=getwt;
Y=X.*W/dot(X,W);

for j=1:length(P)
    for i=1:length(T)
        ckwdot(i,j,:)=ckwyp(P(j)*1013250, T(i), Y, chem);
        

        [wdot(i,j,:)]=getwyp(P(j)*1013250, T(i), Y);
    end
end


for i=1:9
    figure
    surf(P,T,ckwdot(:,:,i)); hold on
    title(sprintf('Species %g',i))
    xlabel('P[atm]')
    ylabel('T[K]')
    zlabel('mol/(cm^3*sec)')
    
    figure
    surf(P,T,wdot(:,:,i)); hold on
    title(sprintf('Species %g',i))
    xlabel('P[atm]')
    ylabel('T[K]')
    zlabel('mol/(cm^3*sec)')
end
