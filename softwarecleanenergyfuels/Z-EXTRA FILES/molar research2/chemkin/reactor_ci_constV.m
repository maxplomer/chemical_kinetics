function reactor_ci_constV
clc;clear;format long e;
tic
% add the path of the toolbox
addpath(strcat(pwd, '/chemkin'));

% load the chemkin workspace
Q=load('GRI3.0\CH4.mat');
chem = Q.chem;
species=Q.species;
RU=Q.RU;%erg/(mol*K)
PA=Q.PA;%1 atm = 1 013 250 erg / (centimeter^3)
W=Q.WT;

% number of variables
KK=length(species);
NV=KK+1;

%system conditions
p0 = 1; % atm
P0 = p0*PA;
T0 = 900;
Phi = 1;

%variables
nY = 1:KK;
nT = KK+1;

% setup the initial conditions - auto calculation of stoich ratio of
% oxdiant to fuel
X0 = zeros(KK,1);
fuel_species={'CH4'};
fuel_ratio=[1];
oxidant_species='O2';
dilutant_species='N2';
dilutant_ratio=3.76; %ratio of dilutant to oxidant

%C H O element composition order
c_order=findelement('c'); h_order = findelement('h'); o_order = findelement('o');
c_tot=0;h_tot=0;o_tot=0;
comp=Q.A';
for j=1:length(fuel_species)
    A=findsp(fuel_species{j}, species);
    c_tot = c_tot + fuel_ratio(j) * comp(A,c_order);
    h_tot = h_tot + fuel_ratio(j) * comp(A,h_order);
	o_tot = o_tot + fuel_ratio(j) * comp(A,o_order);
end
amount_oxidant=h_tot/4 + c_tot - o_tot/2;%calculate stoich amount of oxidant
%initial mol fractions
for j=1:length(fuel_species)
   X0(findsp(fuel_species{j}, species))=fuel_ratio(j); 
end
X0(findsp(oxidant_species, species))=amount_oxidant/Phi;
X0(findsp(dilutant_species, species))=dilutant_ratio*amount_oxidant/Phi;
X0 = X0/sum(X0);

% % setup the initial conditions - simple
% X0 = zeros(KK,1);
% I_H2 = findsp('H2');
% I_O2 = findsp('O2');
% I_N2 = findsp('N2');
% 
% %initial mol fractions
% X0(I_H2) = Phi*2;
% X0(I_O2) = 1;
% X0(I_N2) = 3.76;
% X0 = X0/sum(X0);



% IC vectors  into the main IC vector
y0 = zeros(NV,1);
y0(nY) = X0*P0/(RU*T0); %mol frac ->c_i
y0(nT) = T0; %temperature

%time span seconds
tspan = [0 30]; 

% setup and run the matlab integration solver
options = odeset('RelTol', 1.0E-5, 'AbsTol', 1.0E-9);
[t ,y] = ode15s(@ignfun, tspan, y0, options);



% plot the results
 figure

plot(t, y(:, nT))
xlabel('Time (sec)');
ylabel('Temperature (K)');
xlim(tspan);

%species profiles

% figure
% for i=1:9
% semilogx(t, y(:, i))
% xlabel('Time (sec)');
% ylabel('mol concentration');
% xlim(tspan);
% hold on;
% end

toc



    % RHS functions of the ODEs
    %   input:  t is the time 
    %           y: the vector of indepdent variables 
    %   output: dydt: the vector of LHS of the ODE
    %
    function dydt = ignfun(t, y)%time, mol concentration
        dydt = zeros(NV,1);

        Y = y(nY);
        T = y(nT);
       
        %cv=getcv(T);%vector of cv mol based
        cv=ckcvml(T, chem);%vector of cv mol based
        
        %ctot=sum(Y);  %can use these two lines with alternate formulation of energy eqn
        %cvm=sum(cv.*Y)/ctot;  %
        
        %u=getu(T);%energy mol specific
        u=ckuml(T, chem);%energy mol specific

        %wdot=getwc(T, Y);%molelar production rate
        wdot=ckwc(T, Y, chem);%molelar production rate
        
        % change in mol concentration per time, dY/dt
        dYdt = wdot;
           
        % change in temperature per time, dT/dt
        dTdt  = -dot(u, wdot)/dot(Y, cv);
        dydt(nY) = dYdt;
        dydt(nT) = dTdt;
    end




end