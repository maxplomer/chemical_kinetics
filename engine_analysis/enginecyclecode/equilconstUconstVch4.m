% Software for Clean Energy Fuels - Calculates reaction rates from fuel mechanisms
% Copyright (C) 2014  Max Plomer
 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Contact info for Max Plomer
% email: maxplomer@gmail.com
% cell: 203-945-8606
function [X,P,T] = equilconstUconstVch4(P0, T0)
    
    global Nch4 commonTch4;
    [Nch4,commonTch4]=getthermodatach4;

    A=getcompositionch4';
    RU=83145100;%erg/(mol*K)
    PA=1013250; %1 atm = 1 013 250 erg / (centimeter^3)
    WT=getwtch4;
    KK=getKKch4;
    MM=length(getelementsch4);
    
    N = zeros(KK,1);

    N(findspch4('CH4')) = 1;
    N(findspch4('O2')) = 2;

    a0=A*N;

    N0=zeros(KK+MM+2,1); %2nd to last variable is temperature, last is pressure, these are my guesses
    N0(KK+MM+1)=log(T0);
    N0(KK+MM+2)=P0;
    
    %get u's
    u_0=getuch4(T0);
    
    utot_0=u_0'*N;
    
    X0=N/sum(N);
    
    WM0=dot(X0,WT);
    
    v0=RU*T0 /(P0*PA*WM0);
    
    N = fsolve(@fun,N0);
    N(1:KK)=exp(N(1:KK));
    
    T=exp(N(1+KK+MM));
    
    P=N(2+KK+MM);
    
    X=N(1:KK)/sum(N(1:KK));

    
    function F=fun(Nlog) %Nlog(KK+1:KK+MM) are lamdas O H C N AR

       N(1:KK)=exp(Nlog(1:KK));
       
       T=exp(Nlog(KK+MM+1));
       
       P=Nlog(KK+MM+2);
       
       X=N(1:KK)/sum(N(1:KK));
       
       g=getgch4(T)/(RU*T) + log(P); %log of( P/(1 atm) ), also divided equation by (RU*T) 
       
       F1 = g +     log(X)  +  A'*Nlog(KK+1:KK+MM);
       
       F2 = A*N(1:KK) - a0;
       
       u=getuch4(T);
       
       utot=u'*N(1:KK);
       
       F3=utot/utot_0-1;
       
       WM=dot(X,WT);
       
       v=RU*T /(P*PA*WM);
       
       F4=v/v0-1;
       
       F = [F1 ; F2 ; F3; F4];

    end

end



