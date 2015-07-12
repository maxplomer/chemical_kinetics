function reactor_perfectlystirred_lsode
global tau;
global h0;
global Y0;
global P;
RU=83145100;%erg/(mol*K)
PA=1013250; %1 atm = 1 013 250 erg / (centimeter^3)
W=getwt;

% number of variables
KK=getKK;     %# of species
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


% setup the initial conditions - auto calculation of stoich ratio of
% oxdiant to fuel
Y0 = zeros(KK,1);
fuel_species={'CH4'};
fuel_ratio=[1];
oxidant_species='O2';
dilutant_species='N2';
dilutant_ratio=3.76; %ratio of dilutant to oxidant

%C H O element composition order
c_order=findelement('c'); h_order = findelement('h'); o_order = findelement('o');
c_tot=0;h_tot=0;o_tot=0;
comp=getcomposition;
for j=1:length(fuel_species)
    A=findsp(fuel_species{j});
    c_tot = c_tot + fuel_ratio(j) * comp(A,c_order);
    h_tot = h_tot + fuel_ratio(j) * comp(A,h_order);
	o_tot = o_tot + fuel_ratio(j) * comp(A,o_order);
end
amount_oxidant=h_tot/4 + c_tot - o_tot/2;%calculate stoich amount of oxidant
%initial mass fractions
for j=1:length(fuel_species)
   Y0(findsp(fuel_species{j}))=W(findsp(fuel_species{j}))*fuel_ratio(j); 
end
Y0(findsp(oxidant_species))=W(findsp(oxidant_species))*amount_oxidant/Phi;
Y0(findsp(dilutant_species))=W(findsp(dilutant_species))*dilutant_ratio*amount_oxidant/Phi;
Y0 = Y0/sum(Y0);




% % setup the initial conditions -simple
% Y0 = zeros(KK,1);
% I_CH4 = findsp('CH4');
% I_O2 = findsp('O2');
% I_N2 = findsp('N2');
% 
% %initial mass fractions
% Y0(I_CH4) = W(I_CH4)*Phi*0.5;
% Y0(I_O2) = W(I_O2)*1;
% Y0(I_N2) = W(I_N2)*3.76;
% Y0 = Y0/sum(Y0);

% here is where i put the different IC vectors 
% into the main IC vector
y0 = zeros(NV,1);
y0(nY) = Y0;
y0(nT) = T; %temperature

h0=geth(T)./W;

%time span 0 to .1 seconds
tspan = .1; % .1 seconds

%  run the octave-lsode integration solver
y = lsode ("ignfun", y0, (t = linspace (0, tspan, 200)'));

% plot the results
plot(t, y(:, nT));
xlabel('Time (sec)');
ylabel('Temperature (K)');
xlim([0 tspan]);






    % RHS functions of the ODEs
    %   input:  t is the time 
    %           y: the vector of indepdent variables 
    %   output: dydt: the vector of LHS of the ODE
    %
    function dydt = ignfun(y,t)  % mass faction,time
	global tau;
	global h0;
	global Y0;
	global P;
	KK=getKK;     
	NV = KK+1;    
	nY = 1:KK;
	nT = KK+1;
        W=getwt;
        RU=83145100;%erg/(mol*K)

        dydt = zeros(NV,1);

        Y = y(nY);
        T = y(nT);
        
        %calculate average molecular weight
        WM = 1/sum(Y./W);

        rho = P*WM/(RU*T);  %density from p=roe*R_g*T

        cp=getcp(T)./W;%vector of cv mass based

        h=geth(T)./W;
        
        cpm = dot(cp, Y);%dot product of cp vector and mass fraction vector is mean cp
        
        c=(Y./W)*rho;
        
        wdot=getwc(T, c);
        
        % change in mass fractions per time, dY/dt
        dYdt = (Y0-Y)/tau + wdot.*W/rho;
        
        % change in temperature per time, dT/dt
        dTdt = sum(Y.*(h0-h))/(tau*cpm)-sum(h.*wdot.*W)/(rho*cpm);
        
        dydt(nY) = dYdt;
        dydt(nT) = dTdt;
        
    end
end
