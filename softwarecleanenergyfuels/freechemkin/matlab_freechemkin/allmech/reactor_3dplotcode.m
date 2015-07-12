function  reactor_3dplotcode
T0=300; %initial temperature
Rmax=70; %max compression ratio
Rmin=10; %min compression ratio
Rnum=31; %# of compression ratios want to test
fmax=200; %max frequency
fmin=30; %min frequency
fnum=31; %# of frequencies want to test

R=Rmin:(Rmax-Rmin)/(Rnum-1):Rmax;
f=fmin:(fmax-fmin)/(fnum-1):fmax;

Tmax=zeros(Rnum,fnum);         %maximum temperature
ignition=zeros(Rnum,fnum);     %ignition yes=1/no=0
ignitiondelay=zeros(Rnum,fnum);%time of max temp, no ignition equalls 0
ignitiondelay2=zeros(Rnum,fnum);%time of max temp minux time of half fuel mass fraction, no ignition equalls 0

parTmax=zeros(Rnum*fnum,1);         
parignition=zeros(Rnum*fnum,1);     
parignitiondelay=zeros(Rnum*fnum,1);
parignitiondelay2=zeros(Rnum*fnum,1);
parR=zeros(Rnum*fnum,1);
parf=zeros(Rnum*fnum,1);

%put R and f matrix into parallel compatable vector
for i=1:Rnum*fnum
    x=floor((i-0.0001)/Rnum)+1; %which R are we on?
    y=i-(x-1)*Rnum; %which f are we on?
    parR(i)=R(x);
    parf(i)=f(y);
    
end

matlabpool(4)
parfor i=1:Rnum*fnum
    i
    x=floor((i-0.0001)/Rnum)+1; %which R are we on?
    y=i-(x-1)*Rnum; %which f are we on?

    [ Temp, Yfuel , t ] = reactor_adiabatic_volfunct_3d(parR(i),parf(i),T0);
    
    [ parTmax(i) , I_Tmax ] = max(Temp);
    
    
    %find if ignition
    Yi=Yfuel(1);
    Yf=Yfuel(length(Yfuel));
    if (.5*Yi > Yf) ;     %ignition only if half of initial is greater than final YCH4    
        parignition(i)=1;      %ignition matrix, set to 1 if ignition.
        parignitiondelay(i)=t(I_Tmax);
        for j=1:length(Yfuel)
            if (Yfuel(j) > 0.5*Yi)
                I_halfYfuel = j;
            end
        end
        parignitiondelay2(i)=t(I_Tmax)-t(I_halfYfuel);
    end

end
matlabpool('close');

%put parallel results into matrix
for i=1:Rnum*fnum
    x=floor((i-0.0001)/Rnum)+1; %which R are we on?
    y=i-(x-1)*Rnum; %which f are we on?
    Tmax(x,y) = parTmax(i);
    ignition(x,y)=parignition(i);
    ignitiondelay(x,y)=parignitiondelay(i);
    ignitiondelay2(x,y)=parignitiondelay2(i);
end

figure(1)
mesh(f,R,ignitiondelay)

ylabel('R, compression ratio')
xlabel('f, frequency (Hz)')
zlabel('time of Tmax (sec) (Note: 0 is no ignition)')

figure(2)
mesh(f,R,ignitiondelay2)

ylabel('R, compression ratio')
xlabel('f, frequency (Hz)')
zlabel('Dt from half fuel mass fraction to Tmax (sec) (Note: 0 is no ignition)')

figure(3)
mesh(f,R,ignition)

ylabel('R, compression ratio')
xlabel('f, frequency')
zlabel('0=no ignition, 1=ignition')

figure(4)
mesh(f,R,Tmax)

ylabel('R, compression ratio')
xlabel('f, frequency')
zlabel('Max Temp (K)')

end