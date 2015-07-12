function [results]=run_reactor_adiabatic_volfunct_3d(R,f,T0)

if exist('input.txt','file')==2
    delete('input.txt')	
end

fid = fopen('input.txt','w');

%test parameters
%R = 14.0
%f = 60
%T0 = 500

fprintf(fid,'  %f  %f  %f',R,f,T0)
fclose(fid);


if exist('output.txt','file')==2
    delete('output.txt')	
end


system('./reactor_adiabatic_volfunct_3d>>output.txt');

[t,Temp,Yfuel]=textread('output.txt','%f %f %f');

%test plot
%figure(1)
%plot(t, Temp)
%xlabel('Time (sec)');
%ylabel('Temperature (K)');
%
%figure(2)
%plot(t, Yfuel)
%xlabel('Time (sec)');
%ylabel('Mass Fraction CH4 (g/g)');

    [ Tmax , I_Tmax ] = max(Temp);
    %find if ignition
    Yi=Yfuel(1);
    Yf=Yfuel(length(Yfuel));
    ignition=0;
    ignitiondelay=0;
    ignitiondelay2=0;
    if (.5*Yi > Yf) ;     %ignition only if half of initial is greater than final YCH4    
        ignition=1;      % set to 1 if ignition.
        ignitiondelay=t(I_Tmax);
        for j=1:length(Yfuel)
            if (Yfuel(j) > 0.5*Yi)
                I_halfYfuel = j;
            end
        end
        ignitiondelay2=t(I_Tmax)-t(I_halfYfuel);
    end

results = [ignitiondelay,ignitiondelay2,ignition,Tmax];
end
