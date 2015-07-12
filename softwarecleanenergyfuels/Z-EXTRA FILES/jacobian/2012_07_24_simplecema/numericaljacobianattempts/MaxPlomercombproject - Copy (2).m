%
% Combustion Project for prof T. Lu
%
%
%   by Max Plomer
%
function ign

clc;clear;format long;

chem = ckinit;
[m_names,k_names]=ckname(chem);
m_names
k_names

W = ckwt(chem);
[RU,~,PA] = ckrp(chem);

% number of variables
KK = chem.KK;     %# of species
NV = KK+3;    %# of variables

%system conditions
p0 = 1; % atm
P0 = p0*PA;
T0 = 800;
Phi = 1;


%variables
nY = 1:KK;
nP = KK+1;
nT = KK+2;
nV = KK+3;
I_H2 = ckfindsp('H2',chem);I_H2=I_H2(1);
I_O2 = ckfindsp('O2',chem);I_O2=I_O2(1);
I_N2 = ckfindsp('N2',chem);I_N2=I_N2(1);

% setup the initial conditions
Y0 = zeros(KK,1);

%initial mass fractions
Y0(I_H2) = W(I_H2)*Phi*2;
Y0(I_O2) = W(I_O2)*1;
Y0(I_N2) = W(I_N2)*3.76;
Y0 = Y0/sum(Y0);

%initial average molecular weight
WM0=1/sum(Y0./W);

%initial specific volume
v0 = RU*T0 / (P0*WM0);

% here is where i put the different IC vectors 
% into the main IC vector
y0 = zeros(NV,1);
y0(nY) = Y0;
y0(nP) = P0; %pressure
y0(nT) = T0; %temperature
y0(nV) = v0;

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

figure
semilogx(t, y(:, 1),'.');
xlabel('Time (sec)');
ylabel('Mass Fraction');
xlim(tspan);

%post process
dy=1e-20;%perturb
for i=1:length(t)
    J=[];
    for j=1:(KK+3)
        for k=1:(KK+3)
            y(i,:)
            g=ignfun(0,y(i,:));  
            g=g(j);
            p=zeros(KK+3,1);
            p(k)=dy;
            gp=ignfun(t(i),(y(i,:)+p));  
            gp=gp(j);
            J(j,k)=[gp-g]/dy;
        end
    end
    
    eigs(i,1:KK+3)=eig(J);
end 
figure
loglog(t,eigs(:,1),t,eigs(:,2),t,eigs(:,3))



%g=ignfun


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

        dvdt=0;%

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
