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

function g=getgch4(T)
%gibbs free energy at standard state(p=1atm)
%units: erg/mol
%g=h - Ts
global Nch4 commonTch4;
RU=83145100; %erg/(mol*K)
%g=geth(T) - T*gets(T);

T2=T*T;T3=T2*T;T4=T3*T;T5=T4*T;
logT=T*(1-log(T));
g_over_R=zeros(getKK,1);
for i=1:length(commonTch4)
    if T>commonTch4(i)
        g_over_R(i)=Nch4(1,i)*logT +Nch4(2,i)*-(T2)/2 +Nch4(3,i)*-(T3)/6 +Nch4(4,i)*-(T4)/12 +Nch4(5,i)*-(T5)/20 +Nch4(6,i)-Nch4(7,i)*T;
    else
        g_over_R(i)=Nch4(8,i)*logT +Nch4(9,i)*-(T2)/2 +Nch4(10,i)*-(T3)/6 +Nch4(11,i)*-(T4)/12 +Nch4(12,i)*-(T5)/20 +Nch4(13,i)-Nch4(14,i)*T;
    end
end

g=g_over_R*RU;



end