%
% Combustion Project for prof T. Lu
%
%
%   by Max Plomer
%
function ign
addpath('chemkin')
Q=load('GRI3.0\CH4.mat');

chem = Q.chem;
W = Q.WT;
RU = Q.RU;
PA = Q.PA;
% number of variables
KK = Q.chem.KK;     %# of species
NV = KK+3;    %# of variables

%system conditions
p0 = 1; % atm
P0 = p0*PA;
T0 = 300;
Phi = 0.3;
f=3.4;%rotation frequency,this is a value I will vary, original=100
R=20;%compression ratio, this is a value I will vary, original=50
omega=2*pi*f;%angular velocity
alpha=(R-1)/(R+1);

%variables
nY = 1:53;
nP = 54;
nT = 55;
nV = 56;
I_CH4 = findsp('CH4',Q.species);
I_O2 = findsp('O2',Q.species);
I_N2 = findsp('N2',Q.species);

% setup the initial conditions
Y0 = zeros(KK,1);

%initial mass fractions
Y0(I_CH4) = W(I_CH4)*Phi*0.5;
Y0(I_O2) = W(I_O2)*1;
Y0(I_N2) = W(I_N2)*3.76;
Y0 = Y0/sum(Y0);

%initial average molecular weight
WM0=1/sum(Y0./W);

%initial specific volume
v0 = Q.RU*T0 / (P0*WM0);

% here is where i put the different IC vectors 
% into the main IC vector
y0 = zeros(NV,1);
y0(nY) = Y0;
y0(nP) = P0; %pressure
y0(nT) = T0; %temperature
y0(nV) = v0;

%time span 0 to 10 cycle
tspan = [0 1/f]; % 1 cycles

% setup and run the matlab integration solver
options = odeset('RelTol', 1.0E-6, 'AbsTol', 1.0E-9);
[t ,y] = ode15s(@ignfun, tspan, y0, options);

% plot the results
figure
plot(t, y(:, nT), t, y(:, 14)*20000, '.');
xlabel('Time (sec)');
ylabel('Temperature (K)');
xlim([0 1/f]);

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

        cv=ckcvml(T,chem)./W;%vector of cv mass based

        cvm = dot(cv, Y);%dot product of cv vector and mass fraction vector is mean cv

        u=ckuml(T,chem)./W;%energy mass specific

        dvdt=-v0*omega*alpha*sin(omega*t)/(1+alpha);

        wdot=ckwyp(P, T, Y, chem);%molelar production rate
        
        % change in mass fractions per time, dY/dt
        dYdt = wdot.*W/rho;
           
        % change in temperature per time, dT/dt
        dTdt = -(P*dvdt + sum(u.*wdot.*W)/rho)/cvm;
       
        %Change in Pressure per time, dP/dt 
        dPdt = (T*RU*sum(dYdt./W)+(RU/WM)*dTdt-P*dvdt)*rho;
        
        dydt(nY) = dYdt;
        dydt(nT) = dTdt;
        dydt(nP) = dPdt;
        dydt(nV) = dvdt;
        
    end
end
