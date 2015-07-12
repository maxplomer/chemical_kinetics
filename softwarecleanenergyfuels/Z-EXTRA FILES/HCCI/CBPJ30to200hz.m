clc;
clear;
R=30;
f=50;
i=zeros(31,31);%ignition yes=1/no=0
d=ones(31,31)*-.3e-3;   %time difference btw ignition and max temp, no ignition equalls -1e-3
Temp=zeros(31,31);% Tmax
f=zeros(31,1);%frequency 30-200
R=zeros(31,1);%compression ratio 10-70
T0=[300, 400, 500];%initial temperature from 300 K to 500 K
for p=1:3
TT0=T0(p)

    for x=1:31
        RR=60*(x-1)/30+10
        R(x)=RR;

        for z=1:31
            ff=170*(z-1)/30+30
            f(z)=ff;

            [ T, Yfuel , t ] = CBPJfun( RR , ff, TT0 );

            %find time of max temp
            t1=length(T);
            for v=1:t1
               if (T(v)==max(T));
                q1=v;         %vector position of max temp time
               end
            end
            Temp(x,z)=T(q1);
            %find time of ignition
            t2=length(Yfuel);
            Yi=Yfuel(1);
            Yf=Yfuel(length(Yfuel));

            if (.5*Yi > Yf) ;     %ignition only if half of initial is greater than final YCH4    
                i(x,z)=1;      %ignition matrix, set to 1 if ignition.
                for s=1:t2
                    if (Yfuel(s) > .5*Yi) ;
                    q2=s;         %vector position of ignitiion time, last vector position that was greater than half of initial YCH4
                    end
                end
                d(x,z) = t(q2)-t(q1);     %time of ignition minus time of max temp
            end

            %output
            i(x,z)
            d(x,z)

        end

    end


figure(1+3*(p-1))
mesh(f,R,d)

ylabel('R, compression ratio')
xlabel('f, frequency (Hz)')
zlabel('Dt from Tmax to ignition (sec) (Note: -.3e-3 plane is no ignition)')

figure(2+3*(p-1))
mesh(f,R,i)

ylabel('R, compression ratio')
xlabel('f, frequency')
zlabel('0=no ignition, 1=ignition')

figure(3+3*(p-1))
mesh(f,R,Temp)

ylabel('R, compression ratio')
xlabel('f, frequency')
zlabel('Max Temp (K)')

end