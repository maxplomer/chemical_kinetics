
function reactor_adiabatic_volfunct

global N commonT;
[N,commonT]=getthermodata;

RU=83145100;%erg/(mol*K)
PA=1013250; %1 atm = 1 013 250 erg / (centimeter^3)
W=getwt;

% number of variables
KK = getKK;     %# of species
NV = KK+3;    %# of variables

%system conditions
p0 = 1; % atm
P0 = p0*PA;
T0 = 500;
Phi = 0.3;
f=194.3;%rotation frequency,this is a value I will vary, original=100
R=16;%compression ratio, this is a value I will vary, original=50
omega=2*pi*f;%angular velocity
alpha=(R-1)/(R+1);

%variables
nY = 1:KK;
nP = KK+1;
nT = KK+2;
nV = KK+3;



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

%time span 0 to 10 cycle
tspan = [0 1/f]; % 1 cycles

% setup and run the matlab integration solver
options = odeset('RelTol', 1.0E-6, 'AbsTol', 1.0E-9);
[t ,y] = ode15s(@ignfun, tspan, y0, options);

% plot the results
figure
plot(t, y(:, nT));
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

        cv=getcv(T)./W;%vector of cv mass based

        cvm = dot(cv, Y);%dot product of cv vector and mass fraction vector is mean cv

        u=getu(T)./W;%energy mass specific

        dvdt=-v0*omega*alpha*sin(omega*t)/(1+alpha);

        c=(Y./W)*rho;
        
        wdot=getwc(T, c);
        
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
