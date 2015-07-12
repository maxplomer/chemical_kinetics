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

function cp=getcp(T)
%cp=cv + RU
%units: erg/(mol*K)     erg=g*cm^2/s^2=1e-7 J
global N commonT;
RU=83145100; %erg/(mol*K)

T2=T*T;T3=T2*T;T4=T3*T;
cp_over_R=zeros(getKK,1);
for i=1:length(commonT)
    if T>commonT(i)
        cp_over_R(i)=N(1,i) +N(2,i)*T +N(3,i)*T2 +N(4,i)*T3 +N(5,i)*T4;
    else
        cp_over_R(i)=N(8,i) +N(9,i)*T +N(10,i)*T2 +N(11,i)*T3 +N(12,i)*T4;
    end
end

cp=cp_over_R*RU;

end