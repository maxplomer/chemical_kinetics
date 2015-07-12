function reactor
% k_names = 
% 
%     'H2'
%     'O2'
%     'H2O'
%     'H'
%     'O'
%     'OH'
%     'HO2'
%     'H2O2'
%     'N2'
clc;clear;format long;
chem=ckinit;
RU=83145100;%erg/(mol*K)
PA=1013250; %1 atm = 1 013 250 erg / (centimeter^3)
W=getwt;
% number of variables
KK=9;
NV=KK+2;

%system conditions
p0 = 1; % atm
P0 = p0*PA;
T0 = 800;
Phi = 1;

%variables
nY = 1:KK;
nP = KK+1;
nT = KK+2;
%nV = KK+3;
I_H2 = 1;
I_O2 = 2;
I_N2 = 9;

% setup the initial conditions
Y0 = zeros(KK,1);

%initial mass fractions
Y0(I_H2) = W(I_H2)*Phi*2;
Y0(I_O2) = W(I_O2)*1;
Y0(I_N2) = W(I_N2)*3.76;
Y0 = Y0/sum(Y0);

%initial average molecular weight
%WM0=1/sum(Y0./W);

%initial specific volume
%v0 = RU*T0 / (P0*WM0);

% here is where i put the different IC vectors 
% into the main IC vector
y0 = zeros(NV,1);
y0(nY) = Y0;
y0(nP) = P0; %pressure
y0(nT) = T0; %temperature
%y0(nV) = v0;

%time span seconds
tspan = [0 20]; 

% setup and run the matlab integration solver
options = odeset('RelTol', 1.0E-6, 'AbsTol', 1.0E-9);
[t ,y] = ode15s(@ignfun, tspan, y0, options);

% plot the results
 figure

semilogx(t, y(:, nT),'.')
xlabel('Time (sec)');
ylabel('Temperature (K)');
xlim(tspan);

% figure
% semilogx(t, y(:, 1),'.');
% xlabel('Time (sec)');
% ylabel('Mass Fraction');
% xlim(tspan);




    % RHS functions of the ODEs
    %   input:  t is the time 
    %           y: the vector of indepdent variables 
    %   output: dydt: the vector of LHS of the ODE
    %
    function dydt = ignfun(t, y)%time, mass faction
        dydt = zeros(NV,1);

        Y = y(nY);
        T = y(nT);
        P = y(nP);

        %calculate average molecular weight
        WM = 1/sum(Y./W);

        rho = P*WM/(RU*T);  %density from p=roe*R_g*T

        cv=getcv(T);%vector of cv mass based
        for i=1:KK
            X(i)=Y(i)*WM/W(i);
            c(i)=X(i)*rho/WM;
        end
        
        u=getu(T);%energy mass specific

        %dvdt=0;

        %wdot=getwyp(P, T, Y);%molelar production rate
        wdot=ckwyp(P, T, Y,chem);%molelar production rate
        
        % change in mass fractions per time, dY/dt
        dYdt = wdot.*W/rho;
           
        % change in temperature per time, dT/dt
        dTdt = -dot(u, wdot)/dot(c, cv);
       
        %Change in Pressure per time, dP/dt 
        dPdt = (T*RU*sum(dYdt./W)+(RU/WM)*dTdt)*rho;
        
        dydt(nY) = dYdt;
        dydt(nT) = dTdt;
        dydt(nP) = dPdt;
        %dydt(nV) = dvdt;
        
    end




end