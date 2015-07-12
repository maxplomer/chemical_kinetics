clc;clear;
chem=ckinit;
P=10*1013250; % erg/cm^3
T=300:100:3000;
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


for i=1:length(T)
    ckwdot(i,:)=ckwyp(P, T(i), Y, chem);
    [wdot(i,:)]=getwyp(P, T(i), Y);
end



for i=1:9
    figure
    plot(T,ckwdot(:,i)); hold on
    plot(T,wdot(:,i),'.'); hold on
    title(sprintf('Species %g',i))
    xlabel('T[K]')
    ylabel('mol/(cm^3*sec)')
end
