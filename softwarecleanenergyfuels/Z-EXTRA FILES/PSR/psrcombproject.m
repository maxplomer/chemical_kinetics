%
% Combustion Project for prof T. Lu
%
%
%   by Max Plomer
%
function ign

addpath(strcat(pwd, '/chemkin'));
Q=load('GRI3.0\CH4.mat');

chem = Q.chem;
W = Q.WT;
RU = Q.RU;
PA = Q.PA;

% number of variables
KK=length(Q.species);     %# of species
NV = KK+1;    %# of variables

%system conditions
p = 1; % atm
P = p*PA;
T = 1200;
Phi = 0.5;
tau=1; %1 second

%variables
nY = 1:53;
nT = 54;

I_CH4 = findsp('CH4',Q.species);
I_O2 = findsp('O2',Q.species);
I_N2 = findsp('N2',Q.species);

%for plot
I_H2 = findsp('H2',Q.species);
I_CO = findsp('CO',Q.species);

% setup the initial conditions
Y0 = zeros(KK,1);

%initial mass fractions
Y0(I_CH4) = W(I_CH4)*Phi*0.5;
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
plot(t, y(:, nT));
xlabel('Time (sec)');
ylabel('Temperature (K)');
xlim([0 .1]);

figure(2)
plot(t, y(:, I_H2))
xlabel('Time (sec)');
ylabel('Mass fraction, H2');
xlim([0 .1]);

figure(3)
plot(t, y(:, I_CO))
xlabel('Time (sec)');
ylabel('Mass fraction, CO');
xlim([0 .1]);

figure(4)
plot(t,y(:, I_CH4))
xlabel('Time (sec)');
ylabel('Mass fraction, CH4');
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
        
        % change in mass fractions per time, dY/dt
        dYdt = (Y0-Y)/tau + wdot.*W/rho;
        
        % change in temperature per time, dT/dt
        dTdt = sum(Y.*(h0-h))/(tau*cpm)-sum(h.*wdot.*W)/(rho*cpm);
        
        dydt(nY) = dYdt;
        dydt(nT) = dTdt;
        
    end
end
