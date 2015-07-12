function  reactor_3dplotcode
hosts = [ "54.187.65.173"; 
    "54.186.99.65"; 
    "54.187.49.105";
    "54.187.132.223";
    "54.187.77.188"];
N=4;%number of slaves
sockets = connect (hosts);


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

%create list of commands
for i=1:Rnum*fnum
    cmds{i} = sprintf('run_reactor_adiabatic_volfunct_3d(%g,%g,%g)',parR(i),parf(i),T0);
end



%get results
i=0;
while i<length(cmds)
    i
    for j=1:N
        if (i+j)>length(cmds)
            break;
        end
        cmd=sprintf('send (%s, sockets(1, :))',cmds{i+j});
        reval (cmd, sockets(j+1, :));
        
    end
    for j=1:N
        if (i+j)>length(cmds)
            break;
        end

        results = recv (sockets(j+1, :));
        parignitiondelay(i+j) = results(1);
        parignitiondelay2(i+j) = results(2);
        parignition(i+j) = results(3);
        parTmax(i+j) = results(4);
    end
    i=i+N;
end

 scloseall (sockets);

%put parallel results into matrix
for i=1:Rnum*fnum
    x=floor((i-0.0001)/Rnum)+1; %which R are we on?
    y=i-(x-1)*Rnum; %which f are we on?
    Tmax(x,y) = parTmax(i);
    ignition(x,y)=parignition(i);
    ignitiondelay(x,y)=parignitiondelay(i);
    ignitiondelay2(x,y)=parignitiondelay2(i);
end

save data.mat

end