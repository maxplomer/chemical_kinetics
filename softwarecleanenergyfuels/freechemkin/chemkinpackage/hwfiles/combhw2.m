% A sample proglram to use the chemkin matlab toolbox and equilibrium
% solver
function [X] = combhw2(P0, T0, Phi)
    
    % add the path of the toolbox
    addpath(strcat(pwd, '/chemkin'));

    
    % load the chemkin workspace
    Q=load('GRI3.0\CH4.mat');
    chem = Q.chem;
    A=Q.A;
    RU=Q.RU;
    KK=length(Q.species);
    MM=length(Q.elements);
    
    N = zeros(KK+MM,1);

    N(findsp('CH4', Q.species)) = Phi/2;
    N(findsp('O2', Q.species)) = 1;
    N(findsp('N2', Q.species)) = 3.76;

    a0=A*N(1:KK);

    N0=zeros(KK+MM+1,1); %last variable is temperature, these are my guesses
    N0(KK+MM+1)=T0;
    %get h's
    h_0=ckhml(T0, chem);
    
    htot_0=h_0'*N(1:KK);
    
    N = fsolve(@fun,N0);
    N(1:KK)=exp(N(1:KK));
    
    N(1+KK+MM);
    
    for i=1:KK
       X(i)=N(i)/sum(N(1:KK));
    end
    
    X=X';
    function F=fun(Nlog) %Nlog(KK+1:KK+MM) are lamdas O H C N AR

       N(1:KK)=exp(Nlog(1:KK));
       
       T=Nlog(KK+MM+1);
       
       g=ckgml(T, chem)/(RU*T) + log(P0); %log of p if p not 1
       %/ru*t0, because i change equation around
       
       F1 = g +     log(N(1:KK)/sum(N(1:KK)))  +  A'*Nlog(KK+1:KK+MM);
       
       F2 = A*N(1:KK) - a0;
       
       h=ckhml(T, chem);
       
       htot=h'*N(1:KK);
       
       F3=htot/htot_0-1;
       
       F = [F1 ; F2 ; F3];
    end

end



