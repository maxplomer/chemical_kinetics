%
% Combustion Project for prof T. Lu
%
%
%   by Max Plomer
%
function ign

clc;clear;format long;
chem = ckinit;
W=getwt;
RU=83145100;%erg/(mol*K)
PA=1013250; %1 atm = 1 013 250 erg / (centimeter^3)

% number of variables
KK=9;     %# of species
NV = KK+1;    %# of variables

%system conditions
p = 1; % atm
P = p*PA;
T = 1200;
Phi = 0.5;
tau=1; %1 second

%variables
nY = 1:KK;
nT = KK+1;

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

% here is where i put the different IC vectors 
% into the main IC vector
y0 = zeros(NV,1);
y0(nY) = Y0;
y0(nT) = T; %temperature

h0=ckhml(T, chem)./W;

%time span 0 to .1 seconds
tspan = [0 .1]; % .1 seconds

% setup and run the matlab integration solver
options = odeset('RelTol', 1.0E-6, 'AbsTol', 1.0E-9);
[t ,y] = ode15s(@ignfun, tspan, y0, options);

% plot the results
figure(1)

semilogx(t, y(:, nT),'x');
xlabel('Time (sec)');
ylabel('Temperature (K)');
xlim([0 .1]);





    % RHS functions of the ODEs
    %   input:  t is the time 
    %           y: the vector of indepdent variables 
    %   output: dydt: the vector of LHS of the ODE
    %
    function dydt = ignfun(t, y)  %time, mass faction
        dydt = zeros(NV,1);

        Y = y(nY);
        T = y(nT);
        
        %calculate average molecular weight
        WM = 1/sum(Y./W);

        rho = P*WM/(RU*T);  %density from p=roe*R_g*T

        cp=ckcpml(T,chem)./W;%vector of cv mass based

        h=ckhml(T, chem)./W;
        
        cpm = dot(cp, Y);%dot product of cp vector and mass fraction vector is mean cp
        
        wdot=ckwyp(P, T, Y, chem);%molelar production rate
        %wdot=getwyp(P, T, Y);%molelar production rate
        % change in mass fractions per time, dY/dt
        dYdt = (Y0-Y)/tau + wdot.*W/rho;
        
        % change in temperature per time, dT/dt
        dTdt = sum(Y.*(h0-h))/(tau*cpm)-sum(h.*wdot.*W)/(rho*cpm);
        
        dydt(nY) = dYdt;
        dydt(nT) = dTdt;
        
    end
end
