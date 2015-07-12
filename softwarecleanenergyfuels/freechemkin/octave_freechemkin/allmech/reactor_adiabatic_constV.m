function reactor_adiabatic_constV

RU=83145100;%erg/(mol*K)
PA=1013250; %1 atm = 1 013 250 erg / (centimeter^3)
W=getwt;

% number of variables
KK = getKK;     %# of species
NV = KK+2;    %# of variables

%system conditions
p0 = 1; % atm
P0 = p0*PA;
T0 = 900;
Phi = 1;


%variables
nY = 1:KK;
nP = KK+1;
nT = KK+2;




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
y0(nP) = P0; %pressure
y0(nT) = T0; %temperature


%time span 
t = 0:(10/400):(10);% 1 cycles
t = t';

% run the octave-dassl integration solver
dydt=zeros(NV,1);
y = dassl ("ignfun", y0, dydt,t);

% plot the results
figure
plot(t, y(:, nT));
xlabel('Time (sec)');
ylabel('Temperature (K)');


    % RHS functions of the ODEs
    %   input:  t is the time 
    %           y: the vector of indepdent variables 
    %   output: dydt: the vector of LHS of the ODE
    %
    function res = ignfun(y,dydt,t)% mass faction, change mass fraction per time,time


        KK = getKK;   
        NV = KK+2;
        nY = 1:KK;
        nP = KK+1;
        nT = KK+2;

        W=getwt;
        RU=83145100;%erg/(mol*K)

        Y = y(nY);
        T = y(nT);
        P = y(nP);

        %calculate average molecular weight
        WM = 1/sum(Y./W);

        rho = P*WM/(RU*T);  %density from p=roe*R_g*T

        cv=getcv(T)./W;%vector of cv mass based

        cvm = dot(cv, Y);%dot product of cv vector and mass fraction vector is mean cv

        u=getu(T)./W;%energy mass specific
        
        c=(Y./W)*rho;
        
        wdot=getwc(T, c);
        
        % change in mass fractions per time, dY/dt
        res(nY) = wdot.*W/rho - dydt(nY);
           
        % change in temperature per time, dT/dt
        res(nT) = dydt(nT) + (sum(u.*wdot.*W)/rho)/cvm;
       
        %Change in Pressure per time, dP/dt 
        res(nP) = dydt(nP) - (T*RU*sum(dydt(nY)./W)+(RU/WM)*dydt(nT))*rho;
    end



end
