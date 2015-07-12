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

function g=getg(T)
%gibbs free energy at standard state(p=1atm)
%units: erg/mol
%g=h - Ts
global N commonT;
RU=83145100; %erg/(mol*K)
%g=geth(T) - T*gets(T);

T2=T*T;T3=T2*T;T4=T3*T;T5=T4*T;
logT=T*(1-log(T));
g_over_R=zeros(getKK,1);
for i=1:length(commonT)
    if T>commonT(i)
        g_over_R(i)=N(1,i)*logT +N(2,i)*-(T2)/2 +N(3,i)*-(T3)/6 +N(4,i)*-(T4)/12 +N(5,i)*-(T5)/20 +N(6,i)-N(7,i)*T;
    else
        g_over_R(i)=N(8,i)*logT +N(9,i)*-(T2)/2 +N(10,i)*-(T3)/6 +N(11,i)*-(T4)/12 +N(12,i)*-(T5)/20 +N(13,i)-N(14,i)*T;
    end
end

g=g_over_R*RU;



end