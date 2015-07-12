clc;clear;
%number of reactions

N=21;
chem=ckinit;
P=1013250; %erg/cm^3
T=300:100:3000;
X = zeros(9,1);
X =X+0.01;
Phi=1;

%initial mole fractions
X(1) = Phi*2;
X(2) = 1;
X(9) = 3.76;
X = X/sum(X);

for i=1:length(T)
    [ckfwdk(:,i), ckrevk(:,i)] = ckkfkr(P, T(i), X, chem);
    
    [fwdk(:,i),revk(:,i)]=getkfkr(P, T(i), X);
end

for i=1:N
    figure
    semilogy(T,ckfwdk(i,:)); hold on
    semilogy(T,fwdk(i,:),'.'); hold on
    semilogy(T,ckrevk(i,:)); hold on
    semilogy(T,revk(i,:),'x')
    title(sprintf('Reaction %g',i))
    xlabel('T[K]')
    ylabel('mol/sec')
end
