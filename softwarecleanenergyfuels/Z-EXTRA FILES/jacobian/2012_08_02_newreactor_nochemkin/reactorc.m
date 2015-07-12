function reactorc
% k_names =       'H2'     'O2'     'H2O'     'H'     'O'
%                'OH'     'HO2'     'H2O2'     'N2'
clc;clear;format long e;
chem=ckinit;
RU=83145100;%erg/(mol*K)
PA=1013250; %1 atm = 1 013 250 erg / (centimeter^3)
W=getwt;
% number of variables
KK=9;
NV=KK+1;

%system conditions
p0 = 1; % atm
P0 = p0*PA;
T0 = 800;
Phi = 1;

%variables
nY = 1:KK;
nT = KK+1;

I_H2 = 1;
I_O2 = 2;
I_N2 = 9;

% setup the initial conditions
X0 = zeros(KK,1);

%initial mol fractions
X0(I_H2) = Phi*2;
X0(I_O2) = 1;
X0(I_N2) = 3.76;
X0 = X0/sum(X0);



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

plot(t, y(:, nT),'.')
xlabel('Time (sec)');
ylabel('Temperature (K)');
xlim(tspan);

% figure
% for i=1:9
% semilogx(t, y(:, i))
% xlabel('Time (sec)');
% ylabel('mol concentration');
% xlim(tspan);
% hold on;
% end


% C=y(2, nY);
% C=C';
% T=y(2, nT);
% %J=ajactvmax(T,C,chem)
% J=getjacobian(T, C)
% eig(J)
% for i=1:length(y(:, nT))
%     C=y(i, nY);
%     C=C';
%     T=y(i, nT);
%     J=ajactvmax(T,C,chem);
%     %J=getjacobian(T, C);
%     
%     maxeig(i)=max(eig(J));
% end
% 
% save




    % RHS functions of the ODEs
    %   input:  t is the time 
    %           y: the vector of indepdent variables 
    %   output: dydt: the vector of LHS of the ODE
    %
    function dydt = ignfun(t, y)%time, mol concentration
        dydt = zeros(NV,1);

        Y = y(nY);
        T = y(nT);
       
        cv=getcv(T);%vector of cv mol based
        ctot=sum(Y);
        %cvm=sum(cv.*Y)/ctot;
        u=getu(T);%energy mol specific

        %wdot=getwc(T, Y);%molelar production rate
        wdot=ckwc(T, Y,chem);%molelar production rate
        % change in mol concentration per time, dY/dt
        dYdt = wdot;
           
        % change in temperature per time, dT/dt
        dTdt  = -dot(u, wdot)/dot(Y, cv);
        dydt(nY) = dYdt;
        dydt(nT) = dTdt;
    end




end