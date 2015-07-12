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
function [Tf,Pf,v0,v2,work,heat]=decompression_heatrelease_ch4(P0)

global Nch4 commonTch4;
[Nch4,commonTch4]=getthermodatach4;

RU=83145100;%erg/(mol*K)
PA=1013250; %1 atm = 1 013 250 erg / (centimeter^3)
W=getwtch4;
T0=300;

%results from methane constant volume adiabatic equil
[X,Pf,Tf] = equilconstUconstVch4(P0, T0);

WMf = dot(X,W);   %gram/mol

v0 = RU*Tf / (Pf*PA*WMf); %cm^3/gram

tspan = [0 1];

options = odeset('RelTol', 1.0E-6, 'AbsTol', 1.0E-9);
[t ,y] = ode15s(@ignfun, tspan, Tf, options);

T2=y(length(y));
v2=RU*T2 / (P0*PA*WMf);

u1=dot(getuch4(Tf),X)/WMf;
u2=dot(getuch4(T2),X)/WMf;
work=(u1-u2)/(1e7); %Joule /gram   (1e7 erg = 1 J)

h2=dot(gethch4(T2),X)/WMf;
h3=dot(gethch4(T0),X)/WMf;

heat=(h2-h3)/(1e7); %Joule /gram   (1e7 erg = 1 J)



    function dTdt = ignfun(t, T)%time, mass faction
        
        dPdt = (P0-Pf)*PA;

        P=Pf*PA+dPdt*t;

        cp = getcpch4(T);

        cpm = dot(X,cp)/WMf;

        v = RU*T / (P*WMf);

        dTdt = v*dPdt/cpm;


    end

end