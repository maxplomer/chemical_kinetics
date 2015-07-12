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

function h=geth(T)
%enthalpy, h=f(T) for ideal gas, h=f(P,T) for real gas
%units: erg/mol      erg=g*cm^2/s^2=1e-7 J
global N commonT;
RU=83145100; %erg/(mol*K)
%h=integral of cp dT from 298 to T    +     h_0,f,298

T2=T*T;T3=T2*T;T4=T3*T;T5=T4*T;
h_over_R=zeros(getKK,1);
for i=1:length(commonT)
    if T>commonT(i)
        h_over_R(i)=N(1,i)*T +N(2,i)*(T2)/2 +N(3,i)*(T3)/3 +N(4,i)*(T4)/4 +N(5,i)*(T5)/5 +N(6,i);
    else
        h_over_R(i)=N(8,i)*T +N(9,i)*(T2)/2 +N(10,i)*(T3)/3 +N(11,i)*(T4)/4 +N(12,i)*(T5)/5 +N(13,i);
    end
end

h=h_over_R*RU;

end