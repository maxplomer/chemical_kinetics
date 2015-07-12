% Compare equil solver with solver from problem2


clear all

% add the path of the toolbox
addpath(strcat(pwd, '/chemkin'));

% load the chemkin workspace
load GRI3.0/CH4.mat

% 
Phi = 0.5:0.05:1.5; % equivalence ratio

XH2 = zeros(length(Phi),1);
XCO = zeros(length(Phi),1);
XH2O = zeros(length(Phi),1);
XCH4 = zeros(length(Phi),1);
XH2_prob2 = zeros(length(Phi),1);
XCO_prob2 = zeros(length(Phi),1);
XH2O_prob2 = zeros(length(Phi),1);
XCH4_prob2 = zeros(length(Phi),1);

P0 = 1; % 1 atm
T0 = 1000; % 298 K
TEST = 2000; % guessed final temperature, not important



for n=1:length(Phi)
    % set up mixtures for given Phi
    X0 = zeros(chem.KK, 1);
    X0(findsp('CH4', species)) = Phi(n)/2;
    X0(findsp('O2', species)) = 1;
    X0(findsp('N2', species)) = 3.76;
    X0 = X0/sum(X0);
    
    % call EQUIL solver, 
    % option 1: constant temperature, constant pressure, 
    
     sol = equil(int32(5), P0, T0, TEST, X0, chem);%5 condition
    
    X_prob2 = combhw2(P0, T0, Phi(n));
     
    
     
    % record results
    XH2(n)  = sol.X(findsp('H2', species));
    XCO(n) = sol.X(findsp('CO', species));
    XH2O(n)  = sol.X(findsp('H2O', species));
    XCH4(n)  = sol.X(findsp('CH4', species));
    
        % record results from X_prob2
    XH2_prob2(n)  = X_prob2(findsp('H2', species));
    XCO_prob2(n) = X_prob2(findsp('CO', species));
    XH2O_prob2(n)  = X_prob2(findsp('H2O', species));
    XCH4_prob2(n)  = X_prob2(findsp('CH4', species));
    
end


figure(1)
plot(Phi, XH2, Phi, XCO, Phi, XH2O, Phi, XCH4);
hold on
plot(Phi, XH2_prob2,'x', Phi, XCO_prob2, 'o', Phi, XH2O_prob2, '.', Phi, XCH4_prob2, 'v');
title('p = 1 atm, T0 = 1000K');
xlabel('Equivalence ratio');
ylabel('Species mole fraction');


figure(2)
semilogy(Phi, XH2, Phi, XCO, Phi, XH2O, Phi, XCH4);
hold on
semilogy(Phi, XH2_prob2,'x', Phi, XCO_prob2, 'o', Phi, XH2O_prob2, '.', Phi, XCH4_prob2, 'v');
title('p = 1 atm, T0 = 1000K');
xlabel('Equivalence ratio');
ylabel('Species mole fraction');

