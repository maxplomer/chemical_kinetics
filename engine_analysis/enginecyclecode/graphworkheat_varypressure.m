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
function graphworkheat_varypressure

   P0=[0.3 0.5 0.7 1 2 3 4 5 6 7 8 9 10];

   for i=1:length(P0)
       [Tfch4(i),Pfch4(i),v0ch4(i),v2ch4(i),workch4(i),heatch4(i)]=decompression_heatrelease_ch4(P0(i));
       [Tfh2(i),Pfh2(i),v0h2(i),v2h2(i),workh2(i),heath2(i)]=decompression_heatrelease_h2(P0(i));
   end

   
   figure(1)
   plot(P0,Tfch4,'-o')
   hold on;
   plot(P0,Tfh2,'-x')
   xlabel('Initial Pressure (atm)')
   ylabel('Max Temperature (K)')
   title('Max Temperature (K) vs Initial Pressure (atm)')
   
   figure(2)
   plot(P0,Pfch4,'-o')
   hold on;
   plot(P0,Pfh2,'-x')
   xlabel('Initial Pressure (atm)')
   ylabel('Max Pressure (atm)')
   title('Max Pressure (atm) vs Initial Pressure (atm)')
   
   figure(3)
   plot(P0,v2ch4./v0ch4,'-o')
   hold on;
   plot(P0,v2h2./v0h2,'-x')
   xlabel('Initial Pressure (atm)')
   ylabel('Expansion Ratio (cm3/cm3)')
   title('Expansion Ratio (cm3/cm3) vs Initial Pressure (atm)')
   
   figure(4)
   plot(P0,v2ch4-v0ch4,'-o')
   hold on;
   plot(P0,v2h2-v0h2,'-x')
   xlabel('Initial Pressure (atm)')
   ylabel('Change in volume during expansion (cm3)')
   title('Change in volume during expansion (cm3) vs Initial Pressure (atm)')
   
   figure(5)
   plot(P0,workch4,'-o')
   hold on;
   plot(P0,workh2,'-x')
   xlabel('Initial Pressure (atm)')
   ylabel('Work (Joule/gram)')
   title('Work (Joule/gram) vs Initial Pressure (atm)')
   
   figure(6)
   plot(P0,heatch4,'-o')
   hold on;
   plot(P0,heath2,'-x')
   xlabel('Initial Pressure (atm)')
   ylabel('Heat Release (Joule/gram)')
   title('Heat Release (Joule/gram) vs Initial Pressure (atm)')

end



